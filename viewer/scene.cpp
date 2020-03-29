//
// Created by janos on 02.02.20.
//

#include "scene.hpp"
#include "draw_callbacks.hpp"
#include "scene_graph_node.hpp"

#include <Magnum/SceneGraph/Drawable.h>
#include <Magnum/SceneGraph/Scene.h>
#include <Magnum/Shaders/VertexColor.h>
#include <Magnum/Shaders/Phong.h>
#include <Magnum/Shaders/MeshVisualizer.h>
#include <Magnum/Shaders/Flat.h>
#include <Magnum/Image.h>
#include <Magnum/ImageView.h>
#include <Magnum/Primitives/Capsule.h>


#include <Corrade/Utility/Algorithms.h>

#include <fmt/core.h>
#include <fmt/color.h>

#include <any>
#include <map>

using namespace Magnum;
using namespace Corrade;

using namespace Magnum::Math::Literals;

struct Scene::Impl{

    Impl();

    Object* addObject(std::string name, Object object);

    Object* getObject(std::string_view name);

    SceneGraphNode* addNode(std::string nodeName, std::string_view objectName, ShaderType shader);

    SceneGraphNode* addNode(std::string, std::string_view, std::string_view, ShaderType);

    void setDrawMode(std::string_view node, ShaderType shader);

    void reset() { m_scene.children().clear(); m_nodes.clear(); m_objects.clear(); }

    Scene3D& root();

    SceneGraph::DrawableGroup3D& drawables();

    void setViewportSize(const Vector2i& size);

    Scene3D m_scene;

    SceneGraph::DrawableGroup3D m_drawableGroup;

    std::map<std::string, Object, std::less<>> m_objects;
    std::map<std::string, SceneGraphNode*, std::less<>> m_nodes;

    std::map<ShaderType, std::unique_ptr<GL::AbstractShaderProgram>> m_shaders;
};

Object* Scene::Impl::addObject(
        std::string name,
        Object object)
{
    auto [it, inserted] = m_objects.emplace(std::move(name), std::move(object));
    return inserted ? std::addressof(it->second) : nullptr;
}

Object* Scene::Impl::getObject(std::string_view name){
    auto it = m_objects.find(name);
    if(it == m_objects.end())
        return nullptr;
    else
        return std::addressof(it->second);
}

SceneGraphNode* Scene::Impl::addNode(std::string nodeName, std::string_view objectName, ShaderType shader){
    auto it = m_objects.find(objectName);
    if(it == m_objects.end())
        return nullptr;
    auto& object = it->second;
    auto node = new SceneGraphNode(&m_scene, &m_drawableGroup); //ownership is taking by parent node
    m_nodes.emplace(std::move(nodeName), node);
    node->setDrawCallback(object, m_shaders[shader].get(), shader);
    return node;
}

SceneGraphNode* Scene::Impl::addNode(
        std::string nodeName,
        std::string_view parentName,
        std::string_view objectName,
        ShaderType shader){
    auto* parent = m_nodes.find(parentName)->second;
    auto& object = m_objects.find(objectName)->second;
    auto node = new SceneGraphNode(parent, &m_drawableGroup); //ownership is taking by parent node
    m_nodes.emplace(std::move(nodeName), node);
    node->setDrawCallback(object, m_shaders[shader].get(), shader);
    return node;
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

Object* Scene::addObject(std::string name, Object object)
{
    return m_impl->addObject(std::move(name), std::move(object));
}


Object* Scene::getObject(std::string_view name)
{
    return m_impl->getObject(name);
}

SceneGraphNode* Scene::addNode(std::string node, std::string_view obj, ShaderType shader)
{
    return m_impl->addNode(std::move(node), obj, shader);
}

SceneGraphNode* Scene::addNode(std::string node, std::string_view par, std::string_view obj, ShaderType shader)
{
    return m_impl->addNode(std::move(node), par, obj, shader);
}

void Scene::setDrawMode(std::string_view node, ShaderType shader)
{
   m_impl->setDrawMode(node, shader);
}

void Scene::reset() { m_impl->reset(); }
