
#include "scene_graph_node.hpp"
#include "drawable_data.hpp"

#include <Magnum/Shaders/VertexColor.h>
#include <Magnum/Shaders/Phong.h>
#include <Magnum/Shaders/MeshVisualizer.h>
#include <Magnum/Shaders/Flat.h>

using namespace Magnum;



SceneGraphNode::SceneGraphNode(Object3D* parent, SceneGraph::DrawableGroup3D* group):
    Object3D(parent),
    SceneGraph::Drawable3D{*this, group}
{
}

void SceneGraphNode::draw(const Matrix4& tf, SceneGraph::Camera3D& camera){
    callback(tf, camera);
}

void SceneGraphNode::setDrawCallback(DrawableData& obj, GL::AbstractShaderProgram* shader, ShaderType type)
{
    switch(type){
        case ShaderType::MeshVisualizer :
            callback = MeshVisualizerCallback(obj, *dynamic_cast<Shaders::MeshVisualizer3D*>(shader));
            break;
        case ShaderType::Flat :
        case ShaderType::FlatTextured :
            callback = FlatCallback(obj, *dynamic_cast<Shaders::Flat3D*>(shader));
            break;
        case ShaderType::VertexColor :
            callback = VertexColorCallback(obj, *dynamic_cast<Shaders::VertexColor3D*>(shader));
            break;
        case ShaderType::Phong:
            callback = PhongCallback(obj, *dynamic_cast<Shaders::Phong*>(shader));
            break;
    }
    object = std::addressof(obj);
}
