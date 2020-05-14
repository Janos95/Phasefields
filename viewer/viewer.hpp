//
// Created by janos on 08.11.19.
//

#pragma once

#include "drawables.hpp"
#include "arc_ball_camera.hpp"
#include "primitives.hpp"
#include "subdivision.hpp"
#include "optimization_context.hpp"
#include "solver.hpp"
#include "functional.hpp"
#include "connectedness_constraint.hpp"
#include "modica_mortola.hpp"
#include "problem.hpp"
#include "types.hpp"

#include <Magnum/Platform/Sdl2Application.h>
#include <Magnum/ImGuiIntegration/Context.h>
#include <Magnum/DebugTools/ResourceManager.h>

#include <tbb/task_group.h>
#include <mutex>
#include <atomic>

using namespace Mg::Math::Literals;

namespace Mg = Magnum;
namespace Cr = Corrade;

struct Viewer;

struct OptimizationCallback{
    std::atomic_bool& optimize;
    solver::Status operator()(solver::IterationSummary const&);
};


struct Viewer: public Mg::Platform::Application {

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
    void drawMeshIO();
    bool saveMesh(std::string const&);
    void drawBrushOptions();
    void drawOptimizationContext();
    void makeExclusiveVisualizer(Functional*);
    void drawGradientMetaData(GradientMetaData&, bool&, bool&);
    void drawConnectednessConstraintOptions(ConnectednessMetaData<Mg::Double>&, bool&, bool&);
    void drawShaderOptions();
    void startOptimization();
    void stopOptimization();
    Cr::Containers::Pointer<Functional> makeFunctional(FunctionalType);
    void updateFunctionals(Cr::Containers::Array<Cr::Containers::Pointer<Functional>>&);
    void paint();
    void geodesicSearch();
    void updateInternalDataStructures();
    Mg::Vector3 unproject(Mg::Vector2i const&);

    Mg::ImGuiIntegration::Context imgui{Mg::NoCreate};
    bool trackingMouse = false;

    Mg::UnsignedInt numSubdivisions = 0;

    Cr::Containers::Array<Magnum::Vector3d> vertices;
    Cr::Containers::Array<Magnum::UnsignedInt> indices;

    Mg::Trade::MeshData original{Magnum::MeshPrimitive::Points, 0};
    Mg::Trade::MeshData meshData{Magnum::MeshPrimitive::Points, 0};

    Mg::Trade::MeshData wireframe{Magnum::MeshPrimitive::Points, 0};
    Mg::GL::Mesh wireframeMesh{Mg::NoCreate};
    Object3D * wireframeObject = nullptr;
    PhongDiffuseDrawable * wireframeDrawer = nullptr;

    DrawableType drawableType = DrawableType::PhongDiffuse;
    std::unordered_map<ShaderType, Cr::Containers::Pointer<Mg::GL::AbstractShaderProgram>> shaders;

    ColorMapType colorMapType = ColorMapType::Turbo;
    std::unordered_map<ColorMapType, Magnum::GL::Texture2D> colorMapTextures;

    Magnum::SceneGraph::DrawableGroup3D drawableGroup;
    Scene3D scene;
    Corrade::Containers::Optional<ArcBallCamera> camera;
    Mg::Float near = 0.1f, far = 1000.f;
    Mg::Deg fov = 45._degf;

    Mg::GL::Mesh mesh{Mg::NoCreate};
    Mg::GL::Buffer indexBuffer{Mg::NoCreate}, vertexBuffer{Mg::NoCreate};

    MeshDrawable* drawable = nullptr;
    Object3D* object = nullptr;

    Mg::Color4 clearColor = 0x72909aff_rgbaf;
    std::string expression;

    solver::Problem problem;
    solver::Options options;
    Cr::Containers::Array<Mg::Double> phasefield;
    tbb::task_group g;
    std::mutex mutex;
    std::atomic_bool optimizing = false;

    OptimizationCallback optimizationCallback;
    VisualizationFlags update = {};
    Functional::MetaData* exclusiveVisualizer = nullptr;

    //connectedness vis data
    Cr::Containers::Array<Mg::Color3ub> faceColors;
    Cr::Containers::Pointer<Mg::GL::Texture2D> faceTexture;
    Mg::GL::Texture2D* texture = nullptr;

    SharedRessource<Mg::Double> doubleWellScaling;
    SharedRessource<Mg::Double> dirichletScaling;
    SharedRessource<Mg::Double> connectednessScaling;

    Mg::Double phase = 0;
    Mg::Double targetDist = 0.;
    Mg::Double recursiveFilterFactor = 0.05;
    Mg::Double distStep = 0.1;
    Mg::Double maxDist = 20.f;
    bool brushing = false;
    bool stopPainting = true;
    Cr::Containers::Array<std::pair<double, int>> distances;
    Mg::Vector3d point;

    Mg::DebugTools::ResourceManager manager;
    Drawable* debugDrawable = nullptr;
    bool drawDebug = false;
};