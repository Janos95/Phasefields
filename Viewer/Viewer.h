//
// Created by janos on 08.11.19.
//

#pragma once

#include "ArcBall.h"
#include "primitives.hpp"
#include "Solver.h"
#include "Functional.h"
#include "ModicaMortola.h"
#include "RecursiveProblem.h"
#include "Enums.h"
#include "VisualizationProxy.h"
#include "Tree.h"
#include "Phong.h"
#include "UniqueFunction.h"
#include "FastMarchingMethod.h"
#include "KDTree.h"

#include <Magnum/Platform/Sdl2Application.h>
#include <Magnum/ImGuiIntegration/Context.h>
#include <MagnumPlugins/PrimitiveImporter/PrimitiveImporter.h>
#include <Magnum/Math/Color.h>

namespace Phasefield {

using namespace Mg::Math::Literals;

namespace Mg = Magnum;
namespace Cr = Corrade;

struct Viewer;

struct OptimizationCallback {

    bool optimize = true;

    Solver::Status::Value operator()(Solver::IterationSummary const&);
};


struct Viewer : public Mg::Platform::Application {

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

    void drawMeshIO();

    bool saveMesh(std::string const&);

    bool runOptimization(UniqueFunction<bool()>&& cb);

    void drawBrushOptions();

    void drawOptimizationOptions();

    void drawVisualizationOptions();

    void startOptimization();

    void stopOptimization();

    Functional makeFunctional(FunctionalType::Value);

    void paint();

    void updateInternalDataStructures();

    Mg::Vector3 unproject(Mg::Vector2i const&);

    Mg::ImGuiIntegration::Context imgui{Mg::NoCreate};
    bool trackingMouse = false;

    Mg::UnsignedInt numSubdivisions = 0;

    Cr::Containers::Array<Magnum::Vector3d> vertices;
    Cr::Containers::Array<Magnum::UnsignedInt> indices;

    Mg::Trade::MeshData original{Magnum::MeshPrimitive::Points, 0};
    Mg::Trade::MeshData meshData{Magnum::MeshPrimitive::Points, 0};

    Corrade::Containers::Optional<ArcBall> arcBall;
    Mg::Matrix4 projection;
    Mg::Float near = 0.01f, far = 100.f;
    Mg::Deg fov = 45._degf;

    Mg::GL::Mesh mesh{Mg::NoCreate};
    Mg::GL::Buffer indexBuffer{Mg::NoCreate}, vertexBuffer{Mg::NoCreate};

    Shaders::Phong phong;

    Mg::Color4 clearColor = 0x72909aff_rgbaf;

    Solver::RecursiveProblem problem;
    Solver::Options options;
    Tree tree;
    OptimizationCallback optimizationCallback;
    bool isOptimizing = false;

    //connectedness vis data
    Mg::GL::Texture2D faceTexture{Mg::NoCreate};
    Mg::GL::Texture2D* texture = nullptr;

    SharedRessource<Mg::Double> doubleWellScaling;
    SharedRessource<Mg::Double> dirichletScaling;
    SharedRessource<Mg::Double> connectednessScaling;

    Mg::Double phase = 0;
    Mg::Double targetDist = 0.;
    Mg::Double recursiveFilterFactor = 0.05;
    Mg::Double distStep = 0.01;
    Mg::Double maxDist = 20.f;
    bool brushing = false;
    bool stopPainting = true;
    Cr::Containers::Array<std::pair<double, int>> distances;
    Mg::Vector3d point;

    Mg::GL::Mesh axisMesh{Mg::NoCreate};
    bool drawAxis = false;
    bool drawWireframe = false;

    VisualizationProxy proxy;

    Node* currentNode = nullptr;
    Mg::UnsignedInt colorMapIndex = 0;
    Containers::Array<std::pair<const char*, Magnum::GL::Texture2D>> colorMapTextures;

    Cr::PluginManager::Manager<Mg::Trade::AbstractImporter> manager;
    Mg::Trade::PrimitiveImporter primitiveImporter;

    const int segmentationTag;
    const int phasefieldTag;

    KDTree<Mg::Vector3d> kdtree;
    FastMarchingMethod fastMarchingMethod;
};

}