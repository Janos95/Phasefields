//
// Created by janos on 02.02.20.
//

#pragma once

#include "drawable_data.hpp"
#include "drawables.hpp"

#include <Magnum/SceneGraph/RigidMatrixTransformation3D.h>
#include <Magnum/SceneGraph/Object.h>
#include <Magnum/SceneGraph/Scene.h>
#include <Magnum/Trade/MeshData.h>

#include <memory>
#include <map>

using Object3D = Magnum::SceneGraph::Object<Magnum::SceneGraph::RigidMatrixTransformation3D>;
using Scene3D = Magnum::SceneGraph::Scene<Magnum::SceneGraph::RigidMatrixTransformation3D>;

class Scene {
public:

    Scene();

    Drawable* addNode(std::string nodeName, Object&, DrawableType shader);

    Drawable* setDrawableType(std::string_view nodeName, DrawableType type);

    Object* getNode(std::string_view name);

    [[nodiscard]] bool isDirty() const { return m_dirty; }
    void setClean() { m_dirty = false; }
    void setDirty() { m_dirty = true; }

    Scene3D& root();

    Magnum::SceneGraph::DrawableGroup3D& drawables();

    void setViewportSize(const Magnum::Vector2i& size);

private:
    bool m_dirty = false;

    Scene3D m_scene;

    Magnum::SceneGraph::DrawableGroup3D m_drawableGroup;

    std::map<std::string, Object*, std::less<>> m_nodes;

    std::vector<std::unique_ptr<Magnum::GL::AbstractShaderProgram>> m_shaders;
};
