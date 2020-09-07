//
// Created by janos on 08.11.19.
//

#pragma once

#include "ArcBall.h"
//#include "primitives.hpp"
//#include "Solver.h"
//#include "Functional.h"
//#include "ModicaMortola.h"
//#include "RecursiveProblem.h"
//#include "Enums.h"
//#include "VisualizationProxy.h"
#include "Tree.h"
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

struct OptimizationCallback {

    bool optimize = true;

    //Solver::Status::Value operator()(Solver::IterationSummary const&);
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

    //bool runOptimization(UniqueFunction<bool()>&& cb);

    void drawBrushOptions();

    void drawOptimizationOptions();

    void drawVisualizationOptions();

    void startOptimization();

    void stopOptimization();

    //Functional makeFunctional(FunctionalType::Value);

    void paint();

    void updateInternalDataStructures();

    Mg::Vector3 unproject(Mg::Vector2i const&);

    Mg::ImGuiIntegration::Context imgui{Mg::NoCreate};
    bool trackingMouse = false;

    //UnsignedInt numSubdivisions = 0;

    Mg::GL::Mesh glMesh{Mg::NoCreate};
    Mg::GL::Buffer vertexBuffer{Mg::NoCreate};
    Mg::GL::Buffer indexBuffer{Mg::NoCreate};

    Mesh mesh;
    Mg::Trade::MeshData original{Magnum::MeshPrimitive::Points, 0};
    //Mg::Trade::MeshData meshData{Magnum::MeshPrimitive::Points, 0};

    Optional<ArcBall> arcBall;
    Matrix4 projection;
    Float near = 0.01f, far = 100.f;
    Deg fov = 45._degf;

    //Mg::GL::Mesh mesh{Mg::NoCreate};
    //Mg::GL::Buffer indexBuffer{Mg::NoCreate}, vertexBuffer{Mg::NoCreate};

    Mg::Shaders::Phong phong{Mg::NoCreate};

    Color4 clearColor = 0x72909aff_rgbaf;

    //Solver::RecursiveProblem problem;
    //Solver::Options options;
    Tree tree;
    //OptimizationCallback optimizationCallback;
    //bool isOptimizing = false;

    //connectedness vis data
    //Mg::GL::Texture2D faceTexture{Mg::NoCreate};
    //Mg::GL::Texture2D* texture = nullptr;

    //SharedRessource<Mg::Double> doubleWellScaling;
    //SharedRessource<Mg::Double> dirichletScaling;
    //SharedRessource<Mg::Double> connectednessScaling;

    Double phase = 0;
    Double targetDist = 0.;
    Double recursiveFilterFactor = 0.05;
    Double distStep = 0.01;
    Double maxDist = 20.f;
    bool brushing = false;
    bool stopPainting = true;
    Array<std::pair<double, Vertex>> distances;
    Vector3 point;

    //Mg::GL::Mesh axisMesh{Mg::NoCreate};
    //bool drawAxis = false;
    //bool drawWireframe = false;

    //VisualizationProxy proxy;

    Node* currentNode = nullptr;
    UnsignedInt colorMapIndex = 0;
    Array<std::pair<const char*, Magnum::GL::Texture2D>> colorMapTextures;

    //Cr::PluginManager::Manager<Mg::Trade::AbstractImporter> manager;
    //Mg::Trade::PrimitiveImporter primitiveImporter;

    KDTree<Mg::Vector3> kdtree;
    FastMarchingMethod fastMarchingMethod;
};

}