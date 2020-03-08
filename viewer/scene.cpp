//
// Created by janos on 02.02.20.
//

#include "scene.hpp"
#include "default_callback.hpp"

#include <Magnum/SceneGraph/Drawable.h>
#include <Magnum/SceneGraph/Scene.h>
#include <Magnum/Shaders/VertexColor.h>
#include <Magnum/Shaders/Phong.h>
#include <Magnum/Shaders/MeshVisualizer.h>
#include <Magnum/Shaders/Flat.h>
#include <Magnum/MeshTools/Compile.h>
#include <Magnum/Image.h>
#include <Magnum/ImageView.h>
#include <Magnum/MeshTools/Duplicate.h>
#include <Magnum/MeshTools/GenerateNormals.h>
#include <Magnum/MeshTools/Interleave.h>
#include <Magnum/Primitives/Capsule.h>

#include <Corrade/Utility/Algorithms.h>

#include <variant>
#include <map>

using namespace Magnum;
using namespace Corrade;

using namespace Magnum::Math::Literals;

template< class, class = std::void_t<> >
struct has_viewport_size : std::false_type { };

template< class T >
struct has_viewport_size<
        T,
        /* check if type T has a member setViewportSize that takes a Vector2 */
        std::void_t<decltype(std::declval<T>().setViewportSize(std::declval<Vector2>()))>> : std::true_type { };

template<class T>
constexpr auto has_viewport_size_v = has_viewport_size<T>::value;

struct Scene::Impl{
    Impl();

    bool addObject(std::string, const Trade::MeshData&, CompileFlags flags, const Image2D*);

    Object* getObject(const std::string_view& name);

    Scene3D& root();

    SceneGraph::DrawableGroup3D& drawables();

    void setViewportSize(const Vector2i& size);

    using shader_variant = std::variant<Shaders::Flat3D, Shaders::VertexColor3D, Shaders::MeshVisualizer, Shaders::Phong>;

    Scene3D m_scene;

    SceneGraph::DrawableGroup3D m_drawableGroup;

    std::map<std::string, Object, std::less<>> m_objects;
    std::map<std::string, shader_variant, std::less<>> m_shaders;
};

Scene::Impl::Impl()
{
    m_shaders.emplace("flat_textured", Shaders::Flat3D{Shaders::Flat3D::Flag::Textured});
    m_shaders.emplace("flat", Shaders::Flat3D{});
    m_shaders.emplace("vertex_colored", Shaders::VertexColor3D{});
    m_shaders.emplace("phong", Shaders::Phong{});
    m_shaders.emplace("mesh_vis", Shaders::MeshVisualizer{});
}



Trade::MeshData relayout(const Trade::MeshData& meshData, CompileFlags flags){

    auto ColorA = Trade::MeshAttribute::Color;
    auto PositionA = Trade::MeshAttribute::Position;
    auto NormalA = Trade::MeshAttribute::Normal;

    auto attributes = meshData.attributeData();

    auto generateNormals = meshData.primitive() == MeshPrimitive::Triangles && (flags & (CompileFlag::GenerateFlatNormals|CompileFlag::GenerateSmoothNormals));

    Trade::MeshData generated{MeshPrimitive::Points, 0};
    std::vector<Trade::MeshAttributeData> extra;
    if(generateNormals) {
        CORRADE_ASSERT(meshData.attributeCount(Trade::MeshAttribute::Position),
                       "MeshTools::compile(): the mesh has no positions, can't generate normals", generated);
        /* Right now this could fire only if we have 2D positions, which is
           unlikely; in the future it might fire once packed formats are added */
        CORRADE_ASSERT(meshData.attributeFormat(Trade::MeshAttribute::Position) == VertexFormat::Vector3,
                       "MeshTools::compile(): can't generate normals for" << meshData.attributeFormat(Trade::MeshAttribute::Position) << "positions", generated);

        /* If the data already have a normal array, reuse its location,
           otherwise mix in an extra one */

        if(!meshData.hasAttribute(Trade::MeshAttribute::Normal)) {
            extra.emplace_back(Trade::MeshAttribute::Normal, VertexFormat::Vector3, nullptr);
            /* If we reuse a normal location, expect correct type. Again this won't
               fire now, but might in the future once packed formats are added */
        }
        else{
            CORRADE_ASSERT(meshData.attributeFormat(Trade::MeshAttribute::Normal) == VertexFormat::Vector3,
                           "MeshTools::compile(): can't generate normals into format" << meshData.attributeFormat(Trade::MeshAttribute::Normal), generated);
        }
    }

    if(flags & CompileFlag::AddColorAttribute){
        CORRADE_ASSERT(!meshData.hasAttribute(ColorA), "MeshData already has a color attribute", generated);
        extra.emplace_back(Trade::MeshAttribute::Color, VertexFormat::Vector4, nullptr);
    }

    /* If we want flat normals, we need to first duplicate everything using
       the index buffer. Otherwise just interleave the potential extra
       normal attribute in. */
    if(flags & CompileFlag::GenerateFlatNormals && meshData.isIndexed())
        generated = MeshTools::duplicate(meshData, extra);
    else
        generated = MeshTools::interleave(meshData, extra);

    if(generateNormals) {
        /* Generate the normals. If we don't have the index buffer, we can only
               generate flat ones. */
        if (flags & CompileFlag::GenerateFlatNormals || !meshData.isIndexed())
            MeshTools::generateFlatNormalsInto(
                    generated.attribute<Vector3>(Trade::MeshAttribute::Position),
                    generated.mutableAttribute<Vector3>(Trade::MeshAttribute::Normal));
        else
            MeshTools::generateSmoothNormalsInto(generated.indices(),
                                                 generated.attribute<Vector3>(Trade::MeshAttribute::Position),
                                                 generated.mutableAttribute<Vector3>(Trade::MeshAttribute::Normal));
    }

    return generated;
}

bool Scene::Impl::addObject(
        std::string name,
        const Trade::MeshData& meshdata,
        CompileFlags flags,
        const Image2D* image){

    GL::Buffer indices, vertices;
    if(flags){
        auto generated = relayout(meshdata, flags);
        indices.setData(generated.indexData());
        vertices.setData(generated.vertexData());
    }
    else{
        indices.setData(meshdata.indexData());
        vertices.setData(meshdata.vertexData());
    }
    auto mesh = MeshTools::compile(meshdata, indices, vertices);
    mesh = MeshTools::compile(Primitives::capsule3DSolid(20, 20, 50, 2)); //TODO: remove

    Containers::Optional<GL::Texture2D> texture = Containers::NullOpt;
    if(image)
    {
        texture = GL::Texture2D{};
        texture->setWrapping(GL::SamplerWrapping::ClampToEdge)
                .setMagnificationFilter(GL::SamplerFilter::Linear)
                .setMinificationFilter(GL::SamplerFilter::Linear)
                .setStorage(1, GL::TextureFormat(image->format()), image->size())
                .setSubImage(0, {}, *image);
    }

    Object object{
            std::move(vertices),
            std::move(indices),
            std::move(mesh),
            std::move(texture),
            Color4::cyan(),
            nullptr};

    auto [it, inserted] = m_objects.emplace(std::move(name), std::move(object));
    if(!inserted)
        return false;

    auto& obj = it->second;
    auto& flat = texture ? m_shaders["flat_textured"] : m_shaders["flat"];
    auto cb = std::visit([&](auto& s) -> SceneGraphNode::callback_type {
            if constexpr(std::is_same_v<std::decay_t<decltype(s)>, Shaders::Flat3D>)
                return DefaultCallback(obj, s);
            CORRADE_ASSERT(false, "Impossible code path", {});
        },
        flat);

    obj.node = new SceneGraphNode(&m_scene, cb, &m_drawableGroup); //ownership is taking by parent node
    return true;
}

Object* Scene::Impl::getObject(const std::string_view& name){
    auto it = m_objects.find(name);
    if(it == m_objects.end())
        return nullptr;
    else
        return std::addressof(it->second);
}

Scene3D& Scene::Impl::root(){
    return m_scene;
}

SceneGraph::DrawableGroup3D& Scene::Impl::drawables(){
    return m_drawableGroup;
}

void Scene::Impl::setViewportSize(const Vector2i& size){
    for(auto& [_,shader] : m_shaders) {
        std::visit([size = Vector2(size)](auto &s) {
            if constexpr(has_viewport_size_v < std::remove_reference_t<decltype(s)>>)
                s.setViewportSize(size);
        }, shader);
    }
}

Scene::Scene(): m_impl(std::make_unique<Impl>())
{
}

Scene::~Scene() = default;

bool Scene::addObject(
        std::string name,
        const Trade::MeshData& meshdata,
        CompileFlags flags,
        const Image2D* image)
{
    return m_impl->addObject(std::move(name), meshdata, flags, image);
}

Object* Scene::getObject(const std::string_view& name){
    return m_impl->getObject(name);
}

Scene3D& Scene::root(){
    return m_impl->root();
}

SceneGraph::DrawableGroup3D& Scene::drawables(){
    return m_impl->drawables();
}

void Scene::setViewportSize(const Vector2i& size){
    m_impl->setViewportSize(size);
}


