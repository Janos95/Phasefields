//
// Created by janos on 08.11.19.
//

#pragma once

#include "arc_ball_camera.hpp"
#include "phasefield_data.hpp"
#include "primitives.hpp"
#include "subdivision.hpp"
#include "brush.hpp"
#include "optimization_context.hpp"
#include "shader_options.hpp"

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

using namespace Math::Literals;

using Object3D = Magnum::SceneGraph::Object<Magnum::SceneGraph::MatrixTransformation3D>;
using Scene3D = Magnum::SceneGraph::Scene<Magnum::SceneGraph::MatrixTransformation3D>;

namespace Mg = Magnum;
namespace Cr = Corrade;

struct Viewer: public Magnum::Platform::Application{

    explicit Viewer(int argc, char** argv);

    Mg::GL::Mesh grid{Mg::NoCreate};

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
    void drawShaderOptions();
    void drawPrimitiveOptions();
    void drawBrushOptions();
    void paint();
    void geodesicSearch();
    void finalize();
    Vector3 unproject(Mg::Vector2i const&);

    ImGuiIntegration::Context imgui{NoCreate};
    bool trackingMouse = false;

    int numFaces = 0;
    int numVertices = 0;
    Mg::UnsignedInt numSubdivisions = 0;

    Cr::Containers::Array<Magnum::Vector3d> vertices;
    Cr::Containers::Array<Magnum::UnsignedInt> triangles;
    Cr::Containers::Array<Magnum::Double> phasefield;

    Mg::Trade::MeshData original{Magnum::MeshPrimitive::Points, 0};
    Mg::Trade::MeshData meshData{Magnum::MeshPrimitive::Points, 0};

    ShaderType shader = ShaderType::PhongDiffuse;
    std::unordered_map<ShaderType, Cr::Containers::Pointer<Mg::GL::AbstractShaderProgram>> shaders;

    ColorMapType colorMap = ColorMapType::Turbo;
    std::unordered_map<ColorMapType, Magnum::GL::Texture2D> colorMapTextures;

    Magnum::SceneGraph::DrawableGroup3D drawableGroup;
    Scene3D scene;
    Corrade::Containers::Optional<ArcBallCamera> camera;

    Mg::GL::Mesh mesh{Mg::NoCreate};
    Mg::GL::Buffer indexBuffer{Mg::NoCreate}, vertexBuffer{Mg::NoCreate};

    MeshDrawable* drawable = nullptr;
    Object3D* object = nullptr;
    Object3D *pathManipulator = nullptr;

    Mg::Color4 clearColor = 0x72909aff_rgbaf;
    std::string expression;

    solver::Problem problem;
    solver::Options options;
    std::thread thread;
    std::mutex mutex;
    std::atomic_bool optimize = false;
    std::atomic_bool update = false;

    SharedRessource<Mg::Double> doubleWellScaling;
    SharedRessource<Mg::Double> dirichletScaling;
    SharedRessource<Mg::Double> connectednessScaling;

    Cr::Containers::Array<Mg::Trade::MeshData> paths;
    Mg::Shaders::Flat3D pathShader;
    Mg::Color4 pathColor = Mg::Color4::green();

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