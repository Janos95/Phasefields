
#include "scene_graph_node.hpp"
#include "object.hpp"

using namespace Magnum;



SceneGraphNode::SceneGraphNode(Object3D* parent, callback_type callback, SceneGraph::DrawableGroup3D* group):
    Object3D(parent),
    SceneGraph::Drawable3D{*this, group},
    m_callback(std::move(callback))
{
}

void SceneGraphNode::draw(const Matrix4& tf, SceneGraph::Camera3D& camera){
    m_callback(tf, camera);
}
