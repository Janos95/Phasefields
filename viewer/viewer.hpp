//
// Created by janos on 08.11.19.
//

#pragma once

#include "drawables.hpp"
#include "arc_ball_camera.hpp"
#include "primitives.hpp"
#include "subdivision.hpp"
#include "optimization_context.hpp"
#include "shader_options.hpp"
#include "solver.hpp"
#include "functional.hpp"
#include "problem.hpp"

#include <Corrade/Containers/Reference.h>
#include <Corrade/Containers/Optional.h>
#include <Corrade/Containers/Pointer.h>
#include <Corrade/Containers/GrowableArray.h>

#include <Magnum/Platform/Sdl2Application.h>
#include <Magnum/ImGuiIntegration/Context.h>
#include <Magnum/SceneGraph/MatrixTransformation3D.h>
#include <Magnum/SceneGraph/Object.h>
#include <Magnum/SceneGraph/Scene.h>
#include <Magnum/SceneGraph/Drawable.h>

#include <Magnum/Primitives/Cylinder.h>
#include <Magnum/MeshTools/Compile.h>

#include <thread>
#include <mutex>
#include <atomic>

using namespace Math::Literals;

namespace Mg = Magnum;
namespace Cr = Corrade;

struct Viewer;

template<class T>
struct ConnectednessConstraint;

struct UpdatePhasefield{
    Viewer& viewer;
    solver::Status operator()(solver::IterationSummary const&);
};

struct InstanceData {
    Mg::Matrix4 tf;
    Mg::Matrix3 normalMatrix;
    Mg::Color3 color;
};

struct Paths : Object3D, Drawable {
    explicit Paths(Object3D* parent, SceneGraph::DrawableGroup3D& drawables):
        Object(parent),
        Drawable{*this, &drawables},
        shader{Shaders::Phong::Flag::VertexColor|Shaders::Phong::Flag::InstancedTransformation},
        cylinder{MeshTools::compile(Primitives::cylinderSolid(3,5,0.5f))}
    {
        shader.setAmbientColor(0x111111_rgbf)
              .setSpecularColor(0x330000_rgbf)
              .setLightPosition({10.0f, 15.0f, 5.0f});

        /* cylinder mesh, with an (initially empty) instance buffer */
        cylinder.addVertexBufferInstanced(instanceBuffer, 1, 0,
                                          Shaders::Phong::TransformationMatrix{},
                                          Shaders::Phong::NormalMatrix{},
                                          Shaders::Phong::Color3{});

    }

    void draw(const Matrix4& transformation, SceneGraph::Camera3D& camera) override {
        if(!instanceData) return;
        for(auto& [tf, normalTf, _] : instanceData) {
            tf = transformation * tf;
            normalTf = transformation.normalMatrix() * normalTf;
        }

        instanceBuffer.setData(instanceData, GL::BufferUsage::DynamicDraw);
        cylinder.setInstanceCount(instanceData.size());
        shader.setProjectionMatrix(camera.projectionMatrix())
              .draw(cylinder);
    }

    Shaders::Phong shader;
    GL::Mesh cylinder;
    GL::Buffer instanceBuffer;
    Containers::Array<InstanceData> instanceData;
};

struct Viewer: public Magnum::Platform::Application{

    explicit Viewer(int argc, char** argv);

    void drawEvent() override;
    void tickEvent() override;
    void viewportEvent(ViewportEvent& event) override;
    void keyPressEvent(KeyEvent& event) override;
    void mousePressEvent(MouseEvent& event) override;
    void mouseReleaseEvent(MouseEvent& event) override;
    void mouseMoveEvent(MouseMoveEvent& event) override;
    void mouseScrollEvent(MouseScrollEvent& event) override;
    void keyReleaseEvent(KeyEvent& event) override;
    void textInputEvent(TextInputEvent& event) override;

    void drawSubdivisionOptions();
    void makeDrawableCurrent(DrawableType);
    void drawPrimitiveOptions();
    void drawBrushOptions();
    void drawOptimizationContext();
    bool drawConnectednessConstraintOptions(ConnectednessConstraint<Mg::Double>&, VisualizationFlags&);
    void startOptimization();
    void stopOptimization();
    Cr::Containers::Pointer<Functional> makeFunctional(FunctionalType);
    void updateFunctionals(Cr::Containers::Array<Cr::Containers::Pointer<Functional>>&);
    void paint();
    void geodesicSearch();
    void finalize();
    Mg::GL::Texture2D componentsTexture();
    Mg::GL::Texture2D weightsTexture();
    Vector3 unproject(Mg::Vector2i const&);

    ImGuiIntegration::Context imgui{NoCreate};
    bool trackingMouse = false;

    Mg::UnsignedInt numSubdivisions = 0;

    Cr::Containers::Array<Magnum::Vector3d> vertices;
    Cr::Containers::Array<Magnum::UnsignedInt> indices;

    Mg::Trade::MeshData original{Magnum::MeshPrimitive::Points, 0};
    Mg::Trade::MeshData meshData{Magnum::MeshPrimitive::Points, 0};

    DrawableType drawableType = DrawableType::FaceColored;
    std::unordered_map<ShaderType, Cr::Containers::Pointer<Mg::GL::AbstractShaderProgram>> shaders;

    ColorMapType colorMapType = ColorMapType::Turbo;
    std::unordered_map<ColorMapType, Magnum::GL::Texture2D> colorMapTextures;

    Magnum::SceneGraph::DrawableGroup3D drawableGroup;
    Scene3D scene;
    Corrade::Containers::Optional<ArcBallCamera> camera;

    Mg::GL::Mesh mesh{Mg::NoCreate};
    Mg::GL::Buffer indexBuffer{Mg::NoCreate}, vertexBuffer{Mg::NoCreate};

    MeshDrawable* drawable = nullptr;
    Object3D* object = nullptr;

    Mg::Color4 clearColor = 0x72909aff_rgbaf;
    std::string expression;

    solver::Problem problem;
    solver::Options options;
    std::thread thread;
    std::mutex mutex;
    std::atomic_bool optimizing = false;

    UpdatePhasefield phasefieldCallback;
    VisualizationFlags visFlags = VisualizationFlag::Phasefield;
    Functional* exclusiveVisualizer = nullptr;
    VisualizationFlags update = {};
    Cr::Containers::Array<Magnum::Double> phasefield;
    Cr::Containers::Array<Magnum::Double> parameters;

    Cr::Containers::Array<int> components;
    int numComponents;
    Cr::Containers::Array<Double> ws;
    Mg::GL::Texture2D faceTexture{Mg::NoCreate}; //containing per face colors

    SharedRessource<Mg::Double> doubleWellScaling;
    SharedRessource<Mg::Double> dirichletScaling;
    SharedRessource<Mg::Double> connectednessScaling;

    Cr::Containers::Array<InstanceData> instanceData; //tf from cylinder to path section
    Paths* paths = nullptr;

    Mg::Double phase = 0;
    Mg::Double targetDist;
    Mg::Double recursiveFilterFactor = 0.1;
    Mg::Double distStep = 0.1;
    Mg::Double maxDist = 20.f;
    bool brushing = false;
    bool stopPainting = true;
    Containers::Array<std::pair<double, int>> distances;
    Mg::Vector3d point;
};