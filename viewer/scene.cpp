//
// Created by janos on 02.02.20.
//

#include "scene.hpp"

#include <Magnum/SceneGraph/Drawable.h>
#include <Magnum/Shaders/Phong.h>
#include <Magnum/Shaders/MeshVisualizer.h>
#include <Magnum/Shaders/Flat.h>
#include <Magnum/Primitives/Capsule.h>

#include <Corrade/Utility/Algorithms.h>

#include <map>

using namespace Magnum;
using namespace Corrade;

using namespace Magnum::Math::Literals;


Drawable* Scene::addNode(std::string nodeName, Object& object, DrawableType type){
    if(!object.parent())
        object.setParent(&m_scene);
    Drawable* drawable;
    switch (type) {
        case DrawableType::ColorMapPhong :
            drawable = new ColorMapPhongDrawable(object, *m_shaders[4], &m_drawableGroup);
            break;
        case DrawableType::MeshVisualizer :
            drawable = new MeshVisualizerDrawable(object, *m_shaders[3], &m_drawableGroup);
            break;
        //case NodeType::FlatTextured : node = new ColorMapFlatDrawable(data, *m_shaders[1], parentNode, &m_drawableGroup);
    }
    m_nodes.emplace(std::move(nodeName), &object);
    return drawable;
}

Object* Scene::getNode(std::string_view name){
    auto it = m_nodes.find(name);
    if(it == m_nodes.end())
        return nullptr;
    else
        return it->second;
}

Drawable* Scene::setDrawableType(std::string_view name, DrawableType type){
    auto it = m_nodes.find(name);
    if(it == m_nodes.end())
        return nullptr;
    else{
        auto object = it->second;
        object->features().clear();
        switch (type) {
            case DrawableType::ColorMapPhong :
                new ColorMapPhongDrawable(*object, *m_shaders[4], &m_drawableGroup);
                break;
            case DrawableType::MeshVisualizer :
                new MeshVisualizerDrawable(*object, *m_shaders[3], &m_drawableGroup);
                break;
        }
    }
}

Scene3D& Scene::root(){
    return m_scene;
}

SceneGraph::DrawableGroup3D& Scene::drawables(){
    return m_drawableGroup;
}

void Scene::setViewportSize(const Vector2i& size){
    dynamic_cast<Shaders::MeshVisualizer3D*>(m_shaders[3].get())->setViewportSize(Vector2(size));
}

Scene::Scene() {
    m_shaders.emplace_back(new Shaders::Flat3D{});
    m_shaders.emplace_back(new Shaders::Flat3D{Shaders::Flat3D::Flag::Textured});
    m_shaders.emplace_back(new Shaders::VertexColor3D{});
    m_shaders.emplace_back(new Shaders::MeshVisualizer3D{Shaders::MeshVisualizer3D::Flag::Wireframe | Shaders::MeshVisualizer3D::Flag::NormalDirection});
    m_shaders.emplace_back(new Shaders::Phong{Shaders::Phong::Flag::DiffuseTexture});
}

