//
// Created by janos on 02.02.20.
//

#pragma once

#include "drawable_data.hpp"

#include <Magnum/SceneGraph/RigidMatrixTransformation3D.h>
#include <Magnum/SceneGraph/Object.h>
#include <Magnum/SceneGraph/Scene.h>
#include <Magnum/Trade/MeshData.h>

#include <memory>

using Object3D = Magnum::SceneGraph::Object<Magnum::SceneGraph::RigidMatrixTransformation3D>;
using Scene3D = Magnum::SceneGraph::Scene<Magnum::SceneGraph::RigidMatrixTransformation3D>;

class SceneGraphNode;

enum class ShaderType : Magnum::UnsignedShort {
    Flat = 1,
    FlatTextured = 2,
    VertexColor = 3,
    MeshVisualizer = 4,
    Phong = 5,
};

class Scene {
public:

    Scene();
    ~Scene();

    DrawableData* addObject(std::string name, std::unique_ptr<DrawableData> object);

    DrawableData* getObject(std::string_view name);

    SceneGraphNode* addNode(std::string nodeName, SceneGraphNode*, DrawableData&, ShaderType shader);

    SceneGraphNode* getNode(std::string_view name);

    void setDrawMode(std::string_view node, ShaderType shader);

    void reset();

    [[nodiscard]] bool isDirty() const { return m_dirty; }
    void setClean() { m_dirty = false; }
    void setDirty() { m_dirty = true; }

    Scene3D& root();

    Magnum::SceneGraph::DrawableGroup3D& drawables();

    void setViewportSize(const Magnum::Vector2i& size);

private:
    bool m_dirty;

    //pimpl away map
    struct Impl;
    std::unique_ptr<Impl> m_impl;
};
