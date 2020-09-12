//
// Created by janos on 08.11.19.
//

#pragma once

#include "ArcBall.h"
//#include "primitives.hpp"
#include "Solver.h"
#include "Functional.h"
#include "ModicaMortola.h"
#include "RecursiveProblem.h"
//#include "Enums.h"
#include "VisualizationProxy.h"
#include "Tree.h"
#include "Dijkstra.h"
//#include "UniqueFunction.h"
#include "../Mesh/FastMarchingMethod.h"
#include "KDTree.h"
#include "../Mesh/Mesh.h"
#include "Types.h"

#include <Magnum/Platform/Sdl2Application.h>
#include <Magnum/ImGuiIntegration/Context.h>
//#include <MagnumPlugins/PrimitiveImporter/PrimitiveImporter.h>
#include <Magnum/Math/Color.h>
#include <Magnum/Trade/MeshData.h>
#include <Magnum/Shaders/Phong.h>

namespace Phasefield {

using namespace Mg::Math::Literals;

namespace Mg = Magnum;
namespace Cr = Corrade;

struct Viewer;

//struct OptimizationCallback {
//
//    bool optimize = true;
//
//    Solver::Status::Value operator()(Solver::IterationSummary const&);
//};


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

    //void drawSubdivisionOptions();

    void drawMeshIO();

    //bool saveMesh(std::string const&);

    void runOptimization(UniqueFunction<bool()>&& cb);

    void drawBrushOptions();

    void drawOptimizationOptions();

    void drawVisualizationOptions();

    Functional makeFunctional(FunctionalType::Value);

    bool drawFunctionals(Array<Functional>&, size_t& id);

    void startBrushing(Vector3 const&);
    void brush();

    void updateInternalDataStructures();

    Vector3 unproject(Vector2i const&);

    Mg::ImGuiIntegration::Context imgui{Mg::NoCreate};
    bool trackingMouse = false;

    //UnsignedInt numSubdivisions = 0;

    Mesh mesh;
    Mg::Trade::MeshData original{Magnum::MeshPrimitive::Points, 0};

    Optional<ArcBall> arcBall;
    Matrix4 projection;
    Float near = 0.01f, far = 100.f;
    Deg fov = 45._degf;

    Mg::GL::Mesh glMesh{Mg::NoCreate};
    Mg::GL::Buffer indexBuffer{Mg::NoCreate}, vertexBuffer{Mg::NoCreate};

    Phong phongColorMap{Mg::NoCreate};
    Phong phongVertexColors{Mg::NoCreate};

    Color4 clearColor = 0x72909aff_rgbaf;

    Solver::RecursiveProblem problem;
    Solver::Options options;
    Tree tree;
    bool isOptimizing = false;
    bool hierarchicalOptimization = false;

    //connectedness vis data
    //Mg::GL::Texture2D faceTexture{Mg::NoCreate};
    //Mg::GL::Texture2D* texture = nullptr;

    SharedRessource<double> doubleWellScaling;
    SharedRessource<double> dirichletScaling;
    SharedRessource<Mg::Double> connectednessScaling;

    Double phase = 1.;
    Double targetDist = 0.;
    Double recursiveFilterFactor = 0.05;
    Double distStep = 0.01;
    Double maxDist = 20.f;
    bool brushingModeEnabled = false;
    bool brushing = false;
    Array<std::pair<double, Vertex>> distances;
    Vector3 point;

    //Mg::GL::Mesh axisMesh{Mg::NoCreate};
    //bool drawAxis = false;
    //bool drawWireframe = false;

    VisualizationProxy proxy;

    size_t currentNode = 0;
    bool drawSegmentation = true;
    UnsignedInt colorMapIndex = 0;
    Array<std::pair<const char*, Magnum::GL::Texture2D>> colorMapTextures;

    //Cr::PluginManager::Manager<Mg::Trade::AbstractImporter> manager;
    //Mg::Trade::PrimitiveImporter primitiveImporter;

    KDTree<Mg::Vector3> kdtree;
    FastMarchingMethod fastMarchingMethod;
};

}