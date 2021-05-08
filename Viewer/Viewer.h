//
// Created by janos on 08.11.19.
//

#pragma once

#include "ArcBall.h"
#include "Solver.h"
#include "Functional.h"
#include "RecursiveProblem.h"
#include "VisualizationProxy.h"
#include "Tree.h"
#include "FastMarchingMethod.h"
#include "Mesh.h"
#include "Types.h"
#include "PlotCallback.h"
#include "Bvh.h"
#include "LbfgsSolver.h"

#include <Magnum/ImGuiIntegration/Context.h>
#include <Magnum/Math/Color.h>
#include <Magnum/Trade/MeshData.h>
#include <Magnum/Shaders/Phong.h>
#include <Magnum/Shaders/MeshVisualizer.h>

#include <Corrade/Containers/Optional.h>
#include <Corrade/Utility/Resource.h>

#ifdef MAGNUM_TARGET_WEBGL
#include <Magnum/Platform/EmscriptenApplication.h>
#include <emscripten/emscripten.h>
#include <emscripten/html5.h>
#else
#include <Magnum/Platform/GlfwApplication.h>
#endif

#ifdef PHASEFIELD_WITH_VIDEO
#include "VideoSaver.h"
#endif

#include <MagnumPlugins/StanfordImporter/StanfordImporter.h>
#include <MagnumPlugins/StanfordSceneConverter/StanfordSceneConverter.h>

namespace Phasefield {

using namespace Mg::Math::Literals;

namespace Mg = Magnum;
namespace Cr = Corrade;

struct Viewer;
struct MultiGestureEvent;
struct ScrollingBuffer;


struct Viewer : public Mg::Platform::Application {

    explicit Viewer(Arguments const&);

    ~Viewer();

    void drawEvent() override;

    void viewportEvent(ViewportEvent& event) override;

    void keyPressEvent(KeyEvent& event) override;

    void mousePressEvent(MouseEvent& event) override;

    void mouseReleaseEvent(MouseEvent& event) override;

    void mouseMoveEvent(MouseMoveEvent& event) override;

    void mouseScrollEvent(MouseScrollEvent& event) override;

    void keyReleaseEvent(KeyEvent& event) override;

    void textInputEvent(TextInputEvent& event) override;

#ifdef MAGNUM_TARGET_WEBGL
    MouseEvent const& convertToMouseEvent(EmscriptenTouchEvent const&);
    MouseMoveEvent const& convertToMove(EmscriptenTouchEvent const&);
    Int touchStartEvent(EmscriptenTouchEvent const*);
    Int touchMoveEvent(EmscriptenTouchEvent const*);
    Int touchEndEvent(EmscriptenTouchEvent const*);
    Int touchCancelEvent(EmscriptenTouchEvent const*);
#endif

    //void drawSubdivisionOptions();

    void loadScene(const char*, const char*);

    void loadExperiment(const char*, const char*, const char*);

    //bool saveMesh(std::string const&);

    void runOptimization(UniqueFunction<bool()>&& cb);

    void refineLeafNodes(UniqueFunction<bool()>&& cb);

    void setCallbacks(UniqueFunction<bool()>&& cb);

    void drawBrushOptions();

    void drawOptimizationOptions();

    void drawVisualizationOptions();

    void drawIO();

    //void drawMeshEdit();

    Functional makeFunctional(FunctionalType::Value);

    template<class T>
    bool drawFunctionals(Array<T>&, size_t& id);

    void startBrushing(Vector3 const&, Vector3 const&);

    void brush();

    void setScalingFactors();

    Vector3 unproject(Vector2i const&, Float depth);

    Vertex intersectWithPcd(Vector3 const& p, Vector3 const& dir);

    void drawErrorPlot();

    bool saveMesh(const char*);
    void loadConfig(Cr::Utility::ConfigurationGroup const&);
    void saveCurrentConfig(const char*);


    Mg::ImGuiIntegration::Context imgui{Mg::NoCreate};
    bool trackingMouse = false;

    //UnsignedInt numSubdivisions = 0;

    Mesh mesh;
    Mg::Trade::MeshData original{Magnum::MeshPrimitive::Points, 0};

    Optional<ArcBall> arcBall;
    Matrix4 projection;
    Float near = 0.01f, far = 100.f;
    Deg fov = 45._degf;

    GL::Mesh glMesh{Mg::NoCreate};
    GL::Buffer indexBuffer{Mg::NoCreate}, vertexBuffer{Mg::NoCreate};

    Phong phongVertexColors{Mg::NoCreate};
    Phong phongColorMap{Mg::NoCreate};

    MeshVisualizer3D meshVis{Mg::NoCreate};
    bool drawWireFrame = false;

    Color4 clearColor = 0x72909aff_rgbaf;

    Solver::RecursiveProblem problem;
    Solver::Options options;
    Tree tree;
    bool isOptimizing = false;
    bool hierarchicalOptimization = false;
    bool initializeLevel = false;
    Node currentNode{0, &tree};

    //connectedness vis data
    //Mg::GL::Texture2D faceTexture{Mg::NoCreate};
    //Mg::GL::Texture2D* texture = nullptr;

    double epsilon = 0.023;
    double dirichletScaling;
    double connectednessScaling;
    double areaPenaltyScaling;
    double doubleWellScaling;
    double yamabeLambdaScaling;

    double kappa = 2.;
    double totalArea = 0; /* I guess better than uninitialized */

    Double phase = 1.;
    Double targetDist = 0.;
    Double recursiveFilterFactor = 0.05;
    Double distStep = 0.01;
    Double maxDist = 20.f;
    bool brushingModeEnabled = false;
    bool brushing = false;
    Array<std::pair<double, Vertex>> distances;
    Vector3 point;
    Float lastDepth = 1;

    //Mg::GL::Mesh axisMesh{Mg::NoCreate};
    //bool drawAxis = false;
    //bool drawWireframe = false;

    VisualizationProxy proxy;

    struct ColorMap { const char* name; GL::Texture2D texture; StaticArrayView<256, const Vector3ub> colors; };
    UnsignedInt colorMapIndex = 0;
    Array<ColorMap> colorMapData;


    bool animate = false;
    bool recording = false;

#ifdef PHASEFIELD_WITH_VIDEO
    VideoSaver videoSaver;
#endif

    Mg::Trade::StanfordImporter stanfordImporter;
    //Cr::Utility::Resource experiments;

    BVHAdapter bvh;
    size_t lastIntersectionIdx = 0;
    FastMarchingMethod fastMarchingMethod;

    bool paused = false;
    size_t t = 0;
    bool showPlot = false;

    Optional<LbfgsSolver> pollingSolver;

    /*for touch */
    bool isPinching = false;
    bool trackingFinger = false;
    bool trackingFingers = false;
    bool trackingForImGui = false;
    double pinchLength;
    Vector2i previousMouseMovePosition{-1};


    Initialization::Value init = Initialization::NORMAL_CLUSTER;

    bool swapColors= false;
    size_t colorIndexToSwap[2];
    size_t swapIndex = 0;
};


}