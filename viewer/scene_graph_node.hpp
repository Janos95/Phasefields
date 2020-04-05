//author janos meny

#pragma once

#include "scene.hpp"
#include "draw_callbacks.hpp"

#include <Magnum/Shaders/Shaders.h>
#include <Magnum/SceneGraph/Drawable.h>
#include <Magnum/SceneGraph/Object.h>
#include <Magnum/SceneGraph/RigidMatrixTransformation3D.h>

#include <functional>

using Object3D = Magnum::SceneGraph::Object<Magnum::SceneGraph::RigidMatrixTransformation3D>;
using Scene3D = Magnum::SceneGraph::Scene<Magnum::SceneGraph::RigidMatrixTransformation3D>;
using Drawable3D = Magnum::SceneGraph::Drawable3D;

class SceneGraphNode : public Object3D, public Drawable3D {
public:

    using callback_type = std::function<void(const Magnum::Matrix4&, Magnum::SceneGraph::Camera3D&)>;

    explicit SceneGraphNode(Object3D* parent, Magnum::SceneGraph::DrawableGroup3D* group);

    void setDrawCallback(DrawableData& obj, Magnum::GL::AbstractShaderProgram* shader, ShaderType type);

    DrawableData* object;
    callback_type callback;

protected:

    void draw(const Magnum::Matrix4& tf, Magnum::SceneGraph::Camera3D& camera);

};


