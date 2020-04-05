//
// Created by janos on 02.02.20.
//

#include "scene.hpp"
#include "scene_graph_node.hpp"

#include <Magnum/SceneGraph/Drawable.h>
#include <Magnum/Shaders/Phong.h>
#include <Magnum/Shaders/MeshVisualizer.h>
#include <Magnum/Shaders/Flat.h>
#include <Magnum/ImageView.h>
#include <Magnum/Primitives/Capsule.h>

#include <Corrade/Utility/Algorithms.h>

#include <any>
#include <map>

using namespace Magnum;
using namespace Corrade;

using namespace Magnum::Math::Literals;

struct Scene::Impl{

    Impl();

    DrawableData* addObject(std::string name, std::unique_ptr<DrawableData> object);

    DrawableData* getObject(std::string_view name);

    SceneGraphNode* addNode(std::string nodeName, SceneGraphNode*, DrawableData&, ShaderType shader);

    SceneGraphNode* getNode(std::string_view name);

    void setDrawMode(std::string_view node, ShaderType shader);

    void reset() { m_scene.children().clear(); m_nodes.clear(); m_objects.clear(); }

    Scene3D& root();

    SceneGraph::DrawableGroup3D& drawables();

    void setViewportSize(const Vector2i& size);

    Scene3D m_scene;

    SceneGraph::DrawableGroup3D m_drawableGroup;

    std::map<std::string, std::unique_ptr<DrawableData>, std::less<>> m_objects;
    std::map<std::string, SceneGraphNode*, std::less<>> m_nodes;

    std::map<ShaderType, std::unique_ptr<GL::AbstractShaderProgram>> m_shaders;
};

DrawableData* Scene::Impl::addObject(
        std::string name,
        std::unique_ptr<DrawableData> object)
{
    auto [it, inserted] = m_objects.emplace(std::move(name), std::move(object));
    return inserted ? it->second.get() : nullptr;
}

DrawableData* Scene::Impl::getObject(std::string_view name){
    auto it = m_objects.find(name);
    if(it == m_objects.end())
        return nullptr;
    else
        return it->second.get();
}

SceneGraphNode* Scene::Impl::addNode(std::string nodeName, SceneGraphNode* parent, DrawableData& object, ShaderType shader){
    Object3D* parentNode = parent ? static_cast<Object3D*>(parent) : &m_scene;
    auto node = new SceneGraphNode(parentNode, &m_drawableGroup); //ownership is taking by parent node
    m_nodes.emplace(std::move(nodeName), node);
    node->setDrawCallback(object, m_shaders[shader].get(), shader);
    return node;
}

SceneGraphNode* Scene::Impl::getNode(std::string_view name){
    auto it = m_nodes.find(name);
    if(it == m_nodes.end())
        return nullptr;
    else
        return it->second;
}

void Scene::Impl::setDrawMode(std::string_view nodeName, ShaderType shader){
    auto node = m_nodes.find(nodeName)->second;
    node->setDrawCallback(*(node->object), m_shaders[shader].get(), shader);
}

Scene3D& Scene::Impl::root(){
    return m_scene;
}

SceneGraph::DrawableGroup3D& Scene::Impl::drawables(){
    return m_drawableGroup;
}

void Scene::Impl::setViewportSize(const Vector2i& size){
    dynamic_cast<Shaders::MeshVisualizer3D*>(m_shaders[ShaderType::MeshVisualizer].get())->setViewportSize(Vector2(size));
}

Scene::Impl::Impl() {
    m_shaders.try_emplace(ShaderType::Flat, new Shaders::Flat3D{});
    m_shaders.try_emplace(ShaderType::FlatTextured, new Shaders::Flat3D{Shaders::Flat3D::Flag::Textured});
    m_shaders.try_emplace(ShaderType::VertexColor, new Shaders::VertexColor3D{});
    m_shaders.try_emplace(ShaderType::MeshVisualizer, new Shaders::MeshVisualizer3D{Shaders::MeshVisualizer3D::Flag::Wireframe | Shaders::MeshVisualizer3D::Flag::NormalDirection});
    m_shaders.try_emplace(ShaderType::Phong, new Shaders::Phong{Shaders::Phong::Flag::VertexColor});
}

Scene::Scene(): m_impl(std::make_unique<Impl>())
{
}

Scene::~Scene() = default;

Scene3D& Scene::root(){
    return m_impl->root();
}

SceneGraph::DrawableGroup3D& Scene::drawables(){
    return m_impl->drawables();
}

void Scene::setViewportSize(const Vector2i& size){
    m_impl->setViewportSize(size);
}

DrawableData* Scene::addObject(std::string name, std::unique_ptr<DrawableData> object){
    return m_impl->addObject(std::move(name), std::move(object));
}

DrawableData* Scene::getObject(std::string_view name){
    return m_impl->getObject(name);
}

SceneGraphNode* Scene::addNode(std::string node, SceneGraphNode* parent, DrawableData& obj, ShaderType shader){
    return m_impl->addNode(std::move(node), parent, obj, shader);
}

SceneGraphNode* Scene::getNode(std::string_view name){
    return m_impl->getNode(name);
}

void Scene::setDrawMode(std::string_view node, ShaderType shader){
   m_impl->setDrawMode(node, shader);
}

void Scene::reset() { m_impl->reset(); }
