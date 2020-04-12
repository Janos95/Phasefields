//
// Created by janos on 27.02.20.
//

#pragma once

#include "drawable_data.hpp"

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
#include <Magnum/SceneGraph/RigidMatrixTransformation3D.h>

using Drawable = Magnum::SceneGraph::Drawable3D;
using DrawableGroup = Magnum::SceneGraph::DrawableGroup3D;
using Object = Magnum::SceneGraph::Object<Magnum::SceneGraph::RigidMatrixTransformation3D>;

enum class DrawableType : Magnum::UnsignedShort {
    MeshVisualizer = 0,
    ColorMapPhong = 1,
    ColorMapFlat = 2
};

class ColorMapFlatDrawable : public Drawable
{
public:

    explicit ColorMapFlatDrawable(Object& object, Magnum::GL::AbstractShaderProgram& shader, DrawableGroup* group):
            Drawable(object, group),
            m_mesh(dynamic_cast<ColorMapDrawableData&>(object).mesh),
            m_shader(dynamic_cast<Magnum::Shaders::Flat3D&>(shader)),
            m_textures(dynamic_cast<ColorMapDrawableData&>(object).textures)
    {
    }

    void draw(const Magnum::Matrix4& transformationMatrix, Magnum::SceneGraph::Camera3D& camera) override {
        m_shader.bindTexture(m_textures[int(m_type)])
                .setTransformationProjectionMatrix(camera.projectionMatrix() * transformationMatrix)
                .draw(m_mesh);
    }

    void setColorMapping(ColorMapType type) { m_type = type; }

private:

    Magnum::GL::Mesh& m_mesh;
    Magnum::Shaders::Flat3D& m_shader;
    Corrade::Containers::ArrayView<Magnum::GL::Texture2D> m_textures;
    ColorMapType m_type = ColorMapType::Turbo;
};

class VertexColorDrawable : public Drawable
{
public:

    explicit VertexColorDrawable(Object& object, Magnum::GL::AbstractShaderProgram& shader, DrawableGroup* group):
        Drawable(object, group),
        m_mesh(dynamic_cast<DrawableData&>(object).mesh),
        m_shader(dynamic_cast<Magnum::Shaders::VertexColor3D&>(shader))
    {
    }

    void draw(const Magnum::Matrix4& transformationMatrix, Magnum::SceneGraph::Camera3D& camera) override {
        m_shader.setTransformationProjectionMatrix(camera.projectionMatrix() * transformationMatrix)
                .draw(m_mesh);
    }

private:
    Magnum::GL::Mesh& m_mesh;
    Magnum::Shaders::VertexColor3D& m_shader;
};

class ColorMapPhongDrawable : public Drawable
{
public:

    explicit ColorMapPhongDrawable(Object& object, Magnum::GL::AbstractShaderProgram& shader, DrawableGroup* group):
        Drawable(object, group),
        m_mesh(dynamic_cast<ColorMapDrawableData&>(object).mesh),
        m_shader(dynamic_cast<Magnum::Shaders::Phong&>(shader)),
        m_textures(dynamic_cast<ColorMapDrawableData&>(object).textures)
    {
    }

    void draw(const Magnum::Matrix4& transformationMatrix, Magnum::SceneGraph::Camera3D& camera) override {
        m_shader.bindDiffuseTexture(m_textures[int(m_type)])
                .setLightPositions({{10.0f, 10.0f, 10.0f}, {-10.0f, -10.0f, -10.0f}})
                .setTransformationMatrix(transformationMatrix)
                .setNormalMatrix(transformationMatrix.normalMatrix())
                .setProjectionMatrix(camera.projectionMatrix())
                .draw(m_mesh);
    }

    void setColorMapping(ColorMapType type) { m_type = type; }

private:

    Magnum::GL::Mesh& m_mesh;
    Magnum::Shaders::Phong& m_shader;
    Corrade::Containers::ArrayView<Magnum::GL::Texture2D> m_textures;
    ColorMapType m_type = ColorMapType::Turbo;
};

class MeshVisualizerDrawable : public Drawable
{
public:

    explicit MeshVisualizerDrawable(Object& object, Magnum::GL::AbstractShaderProgram& shader, DrawableGroup* group):
        Drawable(object, group),
        m_mesh(dynamic_cast<DrawableData&>(object).mesh),
        m_shader(dynamic_cast<Magnum::Shaders::MeshVisualizer3D&>(shader))
    {
    }

    void draw(const Magnum::Matrix4& transformationMatrix, Magnum::SceneGraph::Camera3D& camera) override {
        m_shader.setColor(m_colorMesh)
                .setWireframeColor(m_colorWireframe)
                .setTransformationMatrix(transformationMatrix)
                .setProjectionMatrix(camera.projectionMatrix())
                .setNormalMatrix(transformationMatrix.normalMatrix())
                .draw(m_mesh);
    }

    auto& shader(){
        return m_shader;
    }

private:

    Magnum::GL::Mesh& m_mesh;
    Magnum::Color4 m_colorMesh = Magnum::Color4(0,0,0,1);
    Magnum::Color4 m_colorWireframe = Magnum::Color4(1.);
    Magnum::Shaders::MeshVisualizer3D& m_shader;
};
