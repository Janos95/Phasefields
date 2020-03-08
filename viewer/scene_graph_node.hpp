//author janos meny

#pragma once


#include <Magnum/Shaders/Shaders.h>
#include <Magnum/GL/Texture.h>
#include <Magnum/SceneGraph/Camera.h>
#include <Magnum/SceneGraph/Drawable.h>
#include <Magnum/Magnum.h>
#include <Magnum/Math/Color.h>
#include <Magnum/Math/Matrix4.h>
#include <Magnum/SceneGraph/Object.h>
#include <Magnum/SceneGraph/MatrixTransformation3D.h>
#include <Magnum/Shaders/Flat.h>
#include <Magnum/GL/Mesh.h>
#include <Magnum/SceneGraph/RigidMatrixTransformation3D.h>
#include <Magnum/GL/AbstractShaderProgram.h>


using Object3D = Magnum::SceneGraph::Object<Magnum::SceneGraph::RigidMatrixTransformation3D>;
using Scene3D = Magnum::SceneGraph::Scene<Magnum::SceneGraph::RigidMatrixTransformation3D>;
using Drawable3D = Magnum::SceneGraph::Drawable3D;

struct Object;


class SceneGraphNode : public Object3D, public Drawable3D {
public:

    using callback_type = std::function<void(const Magnum::Matrix4&, Magnum::SceneGraph::Camera3D&)>;

    explicit SceneGraphNode(Object3D* parent, callback_type callback, Magnum::SceneGraph::DrawableGroup3D* group);

    void setDrawCallback(callback_type callback){ m_callback = std::move(callback); }

protected:

    void draw(const Magnum::Matrix4& tf, Magnum::SceneGraph::Camera3D& camera) ;

    callback_type m_callback;
};


