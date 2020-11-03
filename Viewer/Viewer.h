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
#include "KDTree.h"
#include "Mesh.h"
#include "Types.h"
#include "PlotCallback.h"
#include "Bvh.h"

#include <Magnum/ImGuiIntegration/Context.h>
#include <Magnum/Math/Color.h>
#include <Magnum/Trade/MeshData.h>
#include <Magnum/Shaders/Phong.h>
#include <Magnum/Shaders/MeshVisualizer.h>

#include <Corrade/Containers/Optional.h>
#include <Corrade/Utility/Resource.h>

#ifdef MAGNUM_TARGET_WEBGL
#include <Magnum/Platform/EmscriptenApplication.h>
#else
#include <Magnum/Platform/Sdl2Application.h>
#endif

#ifdef PHASEFIELD_WITH_IO
#include "VideoSaver.h"
#include <MagnumPlugins/AssimpImporter/AssimpImporter.h>
#endif

#include <MagnumPlugins/StanfordImporter/StanfordImporter.h>
#include <MagnumPlugins/StanfordSceneConverter/StanfordSceneConverter.h>

namespace Phasefield {

using namespace Mg::Math::Literals;

namespace Mg = Magnum;
namespace Cr = Corrade;

struct Viewer;

// utility structure for realtime plot
struct ScrollingBuffer {
    size_t maxSize;
    size_t offset;
    Array<Vector2> data;

    ScrollingBuffer() {
        maxSize = 2000;
        offset  = 0;
        arrayReserve(data, maxSize);
    }

    void add(float x, float y) {
        if (data.size() < maxSize)
            arrayAppend(data, InPlaceInit, x, y);
        else {
            data[offset] = Vector2(x,y);
            offset =  (offset + 1) % maxSize;
        }
    }
    void clear() {
        if (data.size() > 0) {
            arrayShrink(data);
            offset  = 0;
        }
    }

    size_t size() const { return data.size(); }

};

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

    //void drawSubdivisionOptions();

    void loadScene(const char*, const char*);

    void loadExperiment(const char*);

    //bool saveMesh(std::string const&);

    void runOptimization(UniqueFunction<bool()>&& cb);

    void refineLeafNodes(UniqueFunction<bool()>&& cb);

    void setCallbacks(UniqueFunction<bool()>&& cb);

    void drawBrushOptions();

    void drawOptimizationOptions();

    void drawVisualizationOptions();

    void drawIO();

    void drawMeshEdit();

    void setAreaConstraint(Node node);

    Functional makeFunctional(FunctionalType::Value);

    bool drawFunctionals(Array<Functional>&, size_t& id);

    void startBrushing(Vector3 const&, Vector3 const&);

    void brush();

    void setScalingFactors();

    Vector3 unproject(Vector2i const&, Float depth);

    Vertex intersectWithPcd(Vector3 const& p, Vector3 const& dir);

    void drawErrorPlot();

    bool saveMesh(const char*);
    bool dumpMesh(const char*);

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
    size_t maximumDepth = 1;
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

#ifdef PHASEFIELD_WITH_IO
    VideoSaver videoSaver;
    Mg::Trade::AssimpImporter assimpImporter;
#endif

    Mg::Trade::StanfordImporter stanfordImporter;
    //Cr::Utility::Resource experiments;

    BVHAdapter bvh;
    size_t lastIntersectionIdx = 0;
    FastMarchingMethod fastMarchingMethod;

    bool paused = false;
    Array<ScrollingBuffer> data;
    Array<bool> show;
    size_t t = 0;
    bool showPlot = false;
};

}