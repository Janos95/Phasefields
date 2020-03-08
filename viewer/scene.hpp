//
// Created by janos on 02.02.20.
//

#pragma once

#include "object.hpp"

#include <Magnum/SceneGraph/RigidMatrixTransformation3D.h>
#include <Magnum/SceneGraph/Object.h>
#include <Magnum/Trade/MeshData.h>

#include <memory>

using Object3D = Magnum::SceneGraph::Object<Magnum::SceneGraph::RigidMatrixTransformation3D>;
using Scene3D = Magnum::SceneGraph::Scene<Magnum::SceneGraph::RigidMatrixTransformation3D>;


enum class CompileFlag: Magnum::UnsignedByte {
    GenerateFlatNormals = 1 << 0,
    GenerateSmoothNormals = 1 << 1,
    AddColorAttribute = 1 << 2,
};

using CompileFlags = Corrade::Containers::EnumSet<CompileFlag>;

CORRADE_ENUMSET_OPERATORS(CompileFlags)

class Scene {
public:

    Scene();
    ~Scene();

    bool addObject(
            std::string name,
            const Magnum::Trade::MeshData& meshdata,
            CompileFlags flags = {},
            const Magnum::Image2D* image = nullptr);

    Object* getObject(const std::string_view& name);

    Scene3D& root();

    Magnum::SceneGraph::DrawableGroup3D& drawables();

    void setViewportSize(const Magnum::Vector2i& size);

private:
    struct Impl;
    std::unique_ptr<Impl> m_impl;
};
