//
// Created by janos on 27.02.20.
//

#pragma once

#include "object.hpp"

#include <Magnum/Math/Matrix4.h>
#include <Magnum/SceneGraph/Camera.h>
#include <Magnum/GL/Mesh.h>
#include <Magnum/Math/Color.h>
#include <Magnum/Shaders/Flat.h>
#include <Magnum/Shaders/VertexColor.h>
#include <Magnum/Shaders/Phong.h>
#include <Magnum/Shaders/MeshVisualizer.h>


class FlatCallback
{
public:

    explicit FlatCallback(Object&, Magnum::Shaders::Flat3D&);

    void operator()(const Magnum::Matrix4& transformationMatrix, Magnum::SceneGraph::Camera3D& camera);

private:

    Magnum::GL::Mesh& m_mesh;
    Magnum::GL::Texture2D* m_texture = nullptr;
    Magnum::Color4 m_color;
    Magnum::Shaders::Flat3D& m_shader;
};

class VertexColorCallback
{
public:

    explicit VertexColorCallback(Object& obj, Magnum::Shaders::VertexColor3D& shader): m_mesh(obj.mesh), m_shader(shader){}

    void operator()(const Magnum::Matrix4& transformationMatrix, Magnum::SceneGraph::Camera3D& camera){
        m_shader.setTransformationProjectionMatrix(camera.projectionMatrix() * transformationMatrix)
                .draw(m_mesh);
    }

private:
    Magnum::GL::Mesh& m_mesh;
    Magnum::Shaders::VertexColor3D& m_shader;
};

class PhongCallback
{
public:

    explicit PhongCallback(Object& obj, Magnum::Shaders::Phong& shader): m_mesh(obj.mesh), m_shader(shader){}

    void operator()(const Magnum::Matrix4& transformationMatrix, Magnum::SceneGraph::Camera3D& camera){
        m_shader.bindTextures(nullptr, m_diffuseTexture, m_specularTexture, nullptr)
                .setLightPosition({5.0f, 5.0f, 7.0f})
                .setTransformationMatrix(transformationMatrix)
                .setNormalMatrix(transformationMatrix.normalMatrix())
                .setProjectionMatrix(camera.projectionMatrix())
                .draw(m_mesh);
    }

private:

    Magnum::GL::Texture2D *m_diffuseTexture, *m_specularTexture;
    Magnum::GL::Mesh& m_mesh;
    Magnum::Shaders::Phong& m_shader;
};

class MeshVisualizerCallback
{
public:

    explicit MeshVisualizerCallback(Object& obj, Magnum::Shaders::MeshVisualizer& shader): m_mesh(obj.mesh), m_shader(shader){}

    void operator()(const Magnum::Matrix4& transformationMatrix, Magnum::SceneGraph::Camera3D& camera){
        m_shader.setColor(m_colorMesh)
                .setWireframeColor(m_colorWireframe)
                .setTransformationProjectionMatrix(camera.projectionMatrix()*transformationMatrix)
                .draw(m_mesh);
    }

private:

    Magnum::GL::Mesh& m_mesh;
    Magnum::Color4 m_colorMesh = Magnum::Color4::blue();
    Magnum::Color4 m_colorWireframe = Magnum::Color4::red();
    Magnum::Shaders::MeshVisualizer& m_shader;
};
