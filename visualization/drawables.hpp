//
// Created by janos on 27.02.20.
//

#pragma once

#include "types.hpp"

#include <Magnum/GL/GL.h>
#include <Magnum/Math/Math.h>
#include <Magnum/Shaders/Shaders.h>
#include <Magnum/SceneGraph/Drawable.h>
#include <Magnum/Math/Color.h>

namespace Cr = Corrade;
namespace Mg = Magnum;

enum class DrawableType : Magnum::Int {
    MeshVisualizer = 0,
    PhongDiffuse = 1,
    FlatTextured = 2,
    FaceColored = 3,
};

struct MeshDrawable : Drawable {
    MeshDrawable(Object3D& obj, Mg::GL::Mesh& m, DrawableGroup* group, Mg::GL::Texture2D* t = nullptr);
    Mg::GL::Mesh& mesh;
    Mg::GL::Texture2D* texture;
};

struct FlatDrawable : MeshDrawable
{

    explicit FlatDrawable(Object3D&, Mg::GL::Mesh&, Mg::GL::AbstractShaderProgram&, DrawableGroup*);

    void draw(const Magnum::Matrix4&, Magnum::SceneGraph::Camera3D&) override;

    Mg::Shaders::Flat3D& shader;
    Mg::Color4 color;
};

struct VertexColorDrawable : MeshDrawable
{
    explicit VertexColorDrawable(Object3D& object, Mg::GL::Mesh& m, Mg::GL::AbstractShaderProgram& shader, DrawableGroup* group);

    void draw(const Magnum::Matrix4&, Mg::SceneGraph::Camera3D&) override;

    Magnum::Shaders::VertexColor3D* shader = nullptr;
};


struct PhongDiffuseDrawable : MeshDrawable
{
public:

    explicit PhongDiffuseDrawable(Object3D&, Mg::GL::Mesh&, Mg::GL::AbstractShaderProgram&, DrawableGroup*);

    void draw(const Mg::Matrix4&, Mg::SceneGraph::Camera3D&) override;

    Mg::Shaders::Phong& shader;
};

struct MeshVisualizerDrawable : MeshDrawable
{

    explicit MeshVisualizerDrawable(Object3D&, Mg::GL::Mesh&, Mg::GL::AbstractShaderProgram&, DrawableGroup*);

    void draw(const Magnum::Matrix4&, Magnum::SceneGraph::Camera3D&) override;

    Mg::Float wireframeWidth = 1.f;
    Mg::Float smoothness = 2.f;
    Mg::Color4 color = Mg::Color4(1);
    Mg::Color4 wireframeColor = Mg::Color4(0,0,0,1);
    Mg::Shaders::MeshVisualizer3D& shader;
};

struct FaceColorDrawable : MeshDrawable
{
    explicit FaceColorDrawable(Object3D&, Mg::GL::Mesh&, Mg::GL::AbstractShaderProgram&, DrawableGroup*);

    void draw(const Magnum::Matrix4&, Magnum::SceneGraph::Camera3D&) override;

    Mg::Shaders::MeshVisualizer3D& shader;
    Mg::Float offset = 0., scale = 1.f;
};
