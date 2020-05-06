//
// Created by janos on 27.02.20.
//

#pragma once

#include "types.hpp"

#include <Magnum/Math/Matrix4.h>
#include <Magnum/SceneGraph/Camera.h>
#include <Magnum/GL/Mesh.h>
#include <Magnum/Math/Color.h>
#include <Magnum/SceneGraph/Drawable.h>
#include <Magnum/Shaders/Flat.h>
#include <Magnum/Shaders/VertexColor.h>
#include <Magnum/Shaders/Phong.h>
#include <Magnum/Shaders/MeshVisualizer.h>
#include <Magnum/SceneGraph/Object.h>
#include <Magnum/GL/AbstractShaderProgram.h>
#include <Magnum/SceneGraph/MatrixTransformation3D.h>
#include <Magnum/GL/Texture.h>
#include <Magnum/GL/DefaultFramebuffer.h>

namespace Cr = Corrade;
namespace Mg = Magnum;


enum class DrawableType : Magnum::Int {
    MeshVisualizer = 0,
    PhongDiffuse = 1,
    FlatTextured = 2,
    FaceColored = 3,
};

struct MeshDrawable : Drawable {
    MeshDrawable(Object3D& obj, Mg::GL::Mesh& m, DrawableGroup* group, Mg::GL::Texture2D* t = nullptr):
        Drawable(obj, group),
        mesh(m),
        texture(t)
    {
    }
    Mg::GL::Mesh& mesh;
    Mg::GL::Texture2D* texture;
};

struct FlatDrawable : MeshDrawable
{

    explicit FlatDrawable(Object3D& object, Mg::GL::Mesh& m, Mg::GL::AbstractShaderProgram& shader, DrawableGroup* group):
            MeshDrawable(object, m, group),
            shader(dynamic_cast<Magnum::Shaders::Flat3D&>(shader))
    {
    }

    void draw(const Magnum::Matrix4& transformationMatrix, Magnum::SceneGraph::Camera3D& camera) override {
        if(texture) shader.bindTexture(*texture);
        else shader.setColor(color);
        shader.setTransformationProjectionMatrix(camera.projectionMatrix() * transformationMatrix)
              .draw(mesh);
    }

    Mg::Shaders::Flat3D& shader;
    Mg::Color4 color;
};

struct VertexColorDrawable : MeshDrawable
{

    explicit VertexColorDrawable(Object3D& object, Mg::GL::Mesh& m, Mg::GL::AbstractShaderProgram& shader, DrawableGroup* group):
        MeshDrawable(object, m, group),
        shader(&dynamic_cast<Mg::Shaders::VertexColor3D&>(shader))
    {
    }

    void draw(const Magnum::Matrix4& transformationMatrix, Mg::SceneGraph::Camera3D& camera) override {
        shader->setTransformationProjectionMatrix(camera.projectionMatrix() * transformationMatrix)
              .draw(mesh);
    }

    Magnum::Shaders::VertexColor3D* shader = nullptr;
};


struct PhongDiffuseDrawable : MeshDrawable
{
public:

    explicit PhongDiffuseDrawable(Object3D& object, Mg::GL::Mesh& m, Mg::GL::AbstractShaderProgram& shader, DrawableGroup* group):
        MeshDrawable(object, m, group),
        shader(dynamic_cast<Magnum::Shaders::Phong&>(shader))
    {
    }

    void draw(const Mg::Matrix4& transformationMatrix, Mg::SceneGraph::Camera3D& camera) override {
        if(!texture) return;
        shader.bindDiffuseTexture(*texture)
              .setLightPosition({10.0f, 10.0f, 10.0f})
              .setTransformationMatrix(transformationMatrix)
              .setNormalMatrix(transformationMatrix.normalMatrix())
              .setProjectionMatrix(camera.projectionMatrix())
              .draw(mesh);
    }

    Mg::Shaders::Phong& shader;
};

struct MeshVisualizerDrawable : MeshDrawable
{

    explicit MeshVisualizerDrawable(Object3D& object, Mg::GL::Mesh& m, Mg::GL::AbstractShaderProgram& s, DrawableGroup* group):
        MeshDrawable(object, m, group),
        shader(dynamic_cast<Mg::Shaders::MeshVisualizer3D&>(s))
    {
    }

    void draw(const Magnum::Matrix4& transformationMatrix, Magnum::SceneGraph::Camera3D& camera) override {

        shader.setViewportSize(Mg::Vector2{Mg::GL::defaultFramebuffer.viewport().size()})
              .setTransformationMatrix(transformationMatrix)
              .setProjectionMatrix(camera.projectionMatrix())
              .setWireframeWidth(wireframeWidth)
              .setWireframeColor(wireframeColor)
              .setSmoothness(smoothness)
              .setColor(color)
              .draw(mesh);
    }

    Mg::Float wireframeWidth = 1.f;
    Mg::Float smoothness = 2.f;
    Mg::Color4 color = Mg::Color4(1);
    Mg::Color4 wireframeColor = Mg::Color4(0,0,0,1);
    Mg::Shaders::MeshVisualizer3D& shader;
};


struct FaceColorDrawable : MeshDrawable
{

    explicit FaceColorDrawable(Object3D& object, Mg::GL::Mesh& m, Mg::GL::AbstractShaderProgram& s, DrawableGroup* group):
            MeshDrawable(object, m, group),
            shader(dynamic_cast<Mg::Shaders::MeshVisualizer3D&>(s))
    {
    }

    void draw(const Magnum::Matrix4& transformationMatrix, Magnum::SceneGraph::Camera3D& camera) override {
        if(!texture) return;
        shader.setViewportSize(Mg::Vector2{Mg::GL::defaultFramebuffer.viewport().size()})
              .setTransformationMatrix(transformationMatrix)
              .setProjectionMatrix(camera.projectionMatrix())
              .setColorMapTransformation(offset, scale)
              .bindColorMapTexture(*texture)
              .draw(mesh);
    }


    Mg::Shaders::MeshVisualizer3D& shader;
    Mg::Float offset = 0., scale = 1.f;
};
