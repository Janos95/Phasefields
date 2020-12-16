//
// Created by janos on 08.11.19.
//

#include "Viewer.h"
#include "Enums.h"
#include "Tag.h"
#include "ModicaMortola.h"
#include "ConnectednessConstraint.h"
#include "DiffuseYamabe.h"
#include "C1Functions.h"

#include <ScopedTimer/ScopedTimer.h>

#include <Corrade/Utility/Algorithms.h>
#include <Corrade/Utility/Configuration.h>
#include <Corrade/Utility/Directory.h>
#include <Corrade/PluginManager/Manager.h>
#include <Corrade/PluginManager/PluginMetadata.h>
#include <Corrade/Containers/GrowableArray.h>
#include <Corrade/Containers/StringView.h>
#include <Corrade/Utility/ConfigurationGroup.h>

#include <Magnum/GL/DefaultFramebuffer.h>
#include <Magnum/GL/Renderer.h>
#include <Magnum/GL/PixelFormat.h>
#include <Magnum/GL/TextureFormat.h>
#include <Magnum/GL/Buffer.h>
#include <Magnum/GL/Mesh.h>

#include <Magnum/Math/Matrix4.h>
#include <Magnum/Math/FunctionsBatch.h>

#include <Magnum/Shaders/Phong.h>

#include <Magnum/Trade/MeshData.h>
#include <Magnum/Image.h>
#include <Magnum/MeshTools/Duplicate.h>
#include <Magnum/ImageView.h>
#include <Magnum/PixelFormat.h>
#include <Magnum/DebugTools/ColorMap.h>
#include <Magnum/MeshTools/Compile.h>
#include <Magnum/MeshTools/RemoveDuplicates.h>
#include <Magnum/ImGuiIntegration/Context.hpp>

#include <implot.h>
#include <sstream>


#ifdef PHASEFIELD_WITH_IO
#include <MagnumPlugins/StanfordSceneConverter/StanfordSceneConverter.h>

#include <vtkCellArray.h>
#include <vtkPoints.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkTriangle.h>
#include <vtkCellData.h>
#endif



static void initializeRessources() {
    CORRADE_RESOURCE_INITIALIZE(Viewer_Rsc)
    CORRADE_RESOURCE_INITIALIZE(Experiments_Rsc)
}

namespace ColorMap = Mg::DebugTools::ColorMap;

namespace Phasefield {

/* for plotting */
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

struct Experiment {
    const char* name;
    const char* meshName;
    const char* treeName;
    const char* confName;
};

static Experiment experiments[] = {
        {"Connectedness Constraint", "capsule_high_res.ply", "capsule_high_res_split.bin", "connectedness.conf"},
        {"Spot Area Penalty", "spot_high_res.ply", nullptr, "spot_area.conf"},
        {"Yamabe Energy Sphere", "sphere.ply", "sphere.bin", "sphere.conf"},
        {"Yamabe Energy Sphere Connected", "sphere.ply", "sphere.bin", "sphere_connected.conf"}
};

using namespace Mg::Math::Literals;
using namespace Cr::Containers::Literals;

bool Viewer::saveMesh(char const* path) {
#ifdef PHASEFIELD_WITH_IO
    Mg::PluginManager::Manager<Mg::Trade::AbstractSceneConverter> manager;
    Mg::Trade::StanfordSceneConverter exporter{manager, "StanfordSceneConverter"};

    Mg::Trade::MeshData md = mesh.meshDataView();

    if(!exporter.convertToFile(path, md)) {
        Mg::Error{} << "Cannot save file to " << path;
        return false;
    }
#endif
    return true;
}

bool Viewer::dumpMesh(char const* path) {
#ifdef PHASEFIELD_WITH_IO
    mesh.requireGaussianCurvature();
    size_t n = mesh.vertexCount();
    size_t m = mesh.faceCount();

    auto points = vtkSmartPointer<vtkPoints>::New();
    for (Vertex v : mesh.vertices()) {
        auto& p = v.position();
        points->InsertNextPoint(p.x(), p.y(), p.z());
    }

    auto triangles = vtkSmartPointer<vtkCellArray>::New();
    for (Face f : mesh.faces()) {
        auto triangle = vtkSmartPointer<vtkTriangle>::New();
        size_t i = 0;
        for(Vertex v : f.vertices()) {
            triangle->GetPointIds()->SetId(i++,v.idx);
        }
        triangles->InsertNextCell(triangle);
    }

    auto produceArray = [](const char* name, size_t size) {
        vtkSmartPointer<vtkDoubleArray> array = vtkSmartPointer<vtkDoubleArray>::New();
        array->SetNumberOfValues(size);
        array->SetName(name);
        return array;
    };

    auto curvature = produceArray("Curvature", n);
    auto phasefield = produceArray("Phasefield", n);
    auto scalar = produceArray("Scalar", n);
    auto nodalAreas = produceArray("Nodal Area", n);

    auto faceCurvature = produceArray("Face Curvature", m);

    for(Vertex v : mesh.vertices()) {
        curvature->SetValue(v.idx, mesh.gaussianCurvature[v]);
        phasefield->SetValue(v.idx, currentNode.phasefield()[v]);
        scalar->SetValue(v.idx, mesh.scalar(v));
        nodalAreas->SetValue(v.idx, mesh.integral[v]);
    }

    for(Face f : mesh.faces()) {
        faceCurvature->SetValue(f.idx, mesh.faceCurvature[f]);
    }

    // Create a polydata object and add the points to it.
    auto polyData = vtkSmartPointer<vtkPolyData>::New();
    polyData->SetPoints(points);
    polyData->SetPolys(triangles);

    auto pointData = polyData->GetPointData();
    pointData->AddArray(curvature);
    pointData->AddArray(scalar);
    pointData->AddArray(nodalAreas);
    pointData->AddArray(phasefield);

    auto cellData = polyData->GetCellData();
    cellData->AddArray(faceCurvature);

    // Write the file
    vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    writer->SetFileName(path);
    writer->SetInputData(polyData);

    writer->Write();
    return true;
#else
    return false;
#endif
}

namespace {


Array<Viewer::ColorMap> makeColorMapTextures() {

    Array<Viewer::ColorMap> textures;
    using L = std::initializer_list<std::pair<char const*, StaticArrayView<256, const Vector3ub>>>;
    for(auto&& [name, colorMap] : L{
            {"Turbo",   ColorMap::turbo()},
            {"Magma",   ColorMap::magma()},
            {"Plasma",  ColorMap::plasma()},
            {"Inferno", ColorMap::inferno()},
            {"Viridis", ColorMap::viridis()},
            {"CoolWarm", ColorMap::coolWarmSmooth()}
    }) {
        const Magnum::Vector2i size{Magnum::Int(colorMap.size()), 1};
        const GL::TextureFormat format =
#ifdef MAGNUM_TARGET_WEBGL
        GL::TextureFormat::RGB;
#else
        GL::TextureFormat::RGB8;
#endif

        Mg::GL::Texture2D texture;
        texture.setMinificationFilter(Magnum::SamplerFilter::Linear)
               .setMagnificationFilter(Magnum::SamplerFilter::Linear)
               .setWrapping(Magnum::SamplerWrapping::ClampToEdge)
               .setStorage(1, format, size)
               .setSubImage(0, {}, Mg::ImageView2D{Magnum::PixelFormat::RGB8Srgb, size, colorMap});
        arrayAppend(textures, InPlaceInit, name, std::move(texture), colorMap);
    }
    return textures;
}

}

//Solver::Status::Value OptimizationCallback::operator()(Solver::IterationSummary const&) {
//    if(!optimize) return Solver::Status::USER_ABORTED;
//    return Solver::Status::CONTINUE;
//}


Viewer::Viewer(Arguments const& arguments) :
        Mg::Platform::Application{arguments, Mg::NoCreate},
        fastMarchingMethod(mesh),
        tree{mesh},
        bvh{mesh},
        proxy(*this),
        problem(tree)
        //experiments("experiments-data")
{
    initializeRessources();

#ifdef MAGNUM_TARGET_WEBGL
    {
        constexpr static auto cbStart = [](int, const EmscriptenTouchEvent * event, void* userData) -> Int {
            return static_cast<Viewer*>(userData)->touchStartEvent(event);
        };
        constexpr static auto cbMove = [](int, const EmscriptenTouchEvent * event, void* userData) -> Int {
            return static_cast<Viewer*>(userData)->touchMoveEvent(event);
        };
        constexpr static auto cbEnd = [](int, const EmscriptenTouchEvent * event, void* userData) -> Int {
            return static_cast<Viewer*>(userData)->touchEndEvent(event);
        };
        constexpr static auto cbCancel = [](int, const EmscriptenTouchEvent * event, void* userData) -> Int {
            return static_cast<Viewer*>(userData)->touchCancelEvent(event);
        };
        emscripten_set_touchstart_callback("#canvas", this, false, cbStart);
        emscripten_set_touchmove_callback("#canvas", this, false, cbMove);
        emscripten_set_touchend_callback("#canvas", this, false, cbEnd);
        emscripten_set_touchcancel_callback("#canvas", this, false, cbCancel);
    }
#endif

    {
        //arrayAppend(options.callbacks, InPlaceInit, optimizationCallback);
        currentNode = Node{0, &tree};
        setScalingFactors();
    }

/* Try 8x MSAA, fall back to zero samples if not possible. Enable only 2x
MSAA if we have enough DPI. */
    {
        const Vector2 dpiScaling = this->dpiScaling({});
        Configuration conf;
        conf.setTitle("Phasefield Viewer")
            .setWindowFlags(Configuration::WindowFlag::Resizable);

#ifndef MAGNUM_TARGET_WEBGL
        conf.setSize({1200, 1200}, dpiScaling);
#endif

        GLConfiguration glConf;
        glConf.setSampleCount(dpiScaling.max() < 2.0f ? 8 : 2);
        if(!tryCreate(conf, glConf))
            create(conf, glConf.setSampleCount(0));
    }

    /* setup shaders, color map textures and mesh*/
    {
        ScopedTimer t{"Setting up opengl stuff (shaders, textures etc)", true};
        //loadMesh("/home/janos/meshes/spot.ply", original);
        //original = Mg::Primitives::grid3DSolid({300,300});
        //original = Mg::Primitives::grid3DSolid({2,2});
        //mesh.setFromData(original);

        //fastMarchingMethod.update();
        //tree.update();

        //kdtree = KDTree{arrayCast<const Vector3>(mesh.positions())};

        //for(Color4& c : mesh.colors())
        //    c = Color4::red();

        glMesh = Mg::GL::Mesh{};
        vertexBuffer = Mg::GL::Buffer{GL::Buffer::TargetHint::Array};
        indexBuffer = Mg::GL::Buffer{GL::Buffer::TargetHint::ElementArray};

        glMesh.setPrimitive(Mg::MeshPrimitive::Triangles)
              .setIndexBuffer(indexBuffer, 0, Mg::MeshIndexType::UnsignedInt)
              .addVertexBuffer(vertexBuffer, 0, Phong::Position{}, Phong::Normal{}, Phong::TextureCoordinates{}, Phong::Color4{});

        phongColorMap = Phong{Phong::Flag::DiffuseTexture, 2};
        phongVertexColors = Phong{Phong::Flag::VertexColor, 2};
        meshVis = Mg::Shaders::MeshVisualizer3D{Mg::Shaders::MeshVisualizer3D::Flag::Wireframe | Mg::Shaders::MeshVisualizer3D::Flag::NoGeometryShader};

        colorMapData = makeColorMapTextures();
    }

    !Debug{} << "Loading experiment from resource";
    /* Load experiment from resource*/
    {
        ScopedTimer loadingTimer{"Loading experiment", true};
        auto& exp = experiments[2];
        loadExperiment(exp.meshName, exp.treeName, exp.confName);
    }

    /* Setup ImGui, ImPlot, load a better font */
    {
        Debug{} << "Creating ImGui context";
        ImGui::CreateContext();
        ImGui::StyleColorsDark();

        ImFontConfig fontConfig;
        fontConfig.FontDataOwnedByAtlas = false;
        const Vector2 size = Vector2{windowSize()}/dpiScaling();
        Cr::Utility::Resource rs{"viewer-data"};
        ArrayView<const char> font = rs.getRaw("SourceSansPro-Regular.ttf");
        ImGui::GetIO().Fonts->AddFontFromMemoryTTF(
                const_cast<char*>(font.data()), Int(font.size()), 16.0f*framebufferSize().x()/size.x(), &fontConfig);

        imgui = Mg::ImGuiIntegration::Context(*ImGui::GetCurrentContext(),
                                                  Vector2{windowSize()}/dpiScaling(), windowSize(), framebufferSize());

        ImPlot::CreateContext(); /* init implot context */

        /* Setup proper blending to be used by ImGui */
        GL::Renderer::setBlendEquation(
                GL::Renderer::BlendEquation::Add, GL::Renderer::BlendEquation::Add);
        GL::Renderer::setBlendFunction(
                GL::Renderer::BlendFunction::SourceAlpha,
                GL::Renderer::BlendFunction::OneMinusSourceAlpha);
    }

    !Debug{} << "Setting up arcball and projection";
    /* Setup arcball and projection matrix */
    {
        const Vector3 eye = Vector3::zAxis(5.0f);
        const Vector3 viewCenter;
        const Vector3 up = Vector3::yAxis();
        arcBall.emplace(eye, viewCenter, up, fov, windowSize());
        arcBall->setLagging(0.85f);

        projection = Matrix4::perspectiveProjection(fov, Vector2{framebufferSize()}.aspectRatio(), 0.01f, 100.0f);
    }

    GL::Renderer::enable(Mg::GL::Renderer::Feature::DepthTest);
    GL::Renderer::enable(Mg::GL::Renderer::Feature::FaceCulling);

    /* Start the timer */
#ifndef CORRADE_TARGET_EMSCRIPTEN
    setSwapInterval(1);
#endif

}

void Viewer::loadConfig(Cr::Utility::ConfigurationGroup const& conf) {
    epsilon = conf.value<double>("epsilon");
    options.solver = Solver::Backend::Value(conf.value<int>("solver"));
    options.max_num_iterations = conf.value<size_t>("max_num_iterations");
    hierarchicalOptimization = conf.value<bool>("hierarchicalOptimization");
    kappa = conf.value<double>("kappa");

    arrayResize(problem.objectives, 0);
    arrayResize(problem.constraints, 0);
    size_t functionalCount = 0;

    for(FunctionalType::Value fType : FunctionalType::range) {
        auto size = conf.groupCount(FunctionalType::to_string(fType));
        for(size_t i = 0; i < size; ++i) {
            auto* group = conf.group(FunctionalType::to_string(fType), i);
            Functional f = makeFunctional(fType);
            f.loadParameters(*group);
            if(group->value<bool>("isObjective")) {
                arrayAppend(problem.objectives, {.f = MOVE(f)});
            } else {
                arrayAppend(problem.constraints, MOVE(f));
            }
            ++functionalCount;
        }
    }

    setScalingFactors();
}


void Viewer::saveCurrentConfig(const char* path) {

    Cr::Utility::Configuration conf;

    conf.setValue("epsilon", epsilon);
    conf.setValue("solver", int(options.solver));
    conf.setValue("hierarchicalOptimization", hierarchicalOptimization);
    conf.setValue("max_num_iterations", options.max_num_iterations);
    conf.setValue("kappa", kappa);

    for(auto& [f, hist, show] : problem.objectives) {
        auto* group = conf.addGroup(FunctionalType::to_string(f.functionalType));
        f.saveParameters(*group);
        group->setValue("isObjective", true);
    }

    for(auto& f : problem.constraints) {
        auto* group = conf.addGroup(FunctionalType::to_string(f.functionalType));
        f.saveParameters(*group);
        group->setValue("isObjective", false);
    }

    conf.save(path);
}

void Viewer::loadScene(const char* path, const char* postfix) {
#ifdef PHASEFIELD_WITH_IO
    bool loadedTree = false;
    Optional<Mg::Trade::MeshData> md;

    if(Cr::Utility::Directory::exists(path)) {

        std::string postfixStr{postfix};
        std::string meshPath = Cr::Utility::Directory::join(path, postfixStr + ".ply");
        std::string treePath = Cr::Utility::Directory::join(path ,postfixStr + ".bin");
        std::string confPath = Cr::Utility::Directory::join(path ,postfixStr + ".conf");

        if(Cr::Utility::Directory::exists(treePath)) {
            FILE* fp = std::fopen(treePath.c_str(), "r");
            fseek(fp, 0, SEEK_END);
            size_t size = std::ftell(fp);
            fseek(fp, 0, SEEK_SET);
            Array<char> data{NoInit, size};
            fread(data.data(), sizeof(char), size, fp);
            fclose(fp);

            tree = Tree::deserialize(data, mesh);
            loadedTree = true;
        }

        if(Cr::Utility::Directory::exists(meshPath)) {

            Debug{} << "Opening file" << meshPath.c_str();

            if(assimpImporter.openFile(meshPath)) {
                Debug{} << "Imported " << assimpImporter.meshCount() << " meshes";

                if(assimpImporter.meshCount() && assimpImporter.mesh(0)) {
                    md = assimpImporter.mesh(0);
                } else {
                    Debug{} << "Could not load mesh";
                }
            } else {
                puts("could not open file");
            }
        }
        if(Cr::Utility::Directory::exists(confPath)) {
            Cr::Utility::Configuration conf{confPath};
            loadConfig(conf);
        }
    }

    if(md) {
        auto cleaned = Mg::MeshTools::removeDuplicatesFuzzy(*md);
        mesh.setFromData(cleaned);

        mesh.uploadVertexBuffer(vertexBuffer);
        mesh.uploadIndexBuffer(indexBuffer);

        glMesh.setPrimitive(Mg::MeshPrimitive::Triangles)
              .setCount(mesh.indexCount())
              .setIndexBuffer(indexBuffer, 0, Mg::MeshIndexType::UnsignedInt)
              .addVertexBuffer(vertexBuffer, 0, Phong::Position{}, Phong::Normal{}, Phong::TextureCoordinates{},
                               Phong::Color4{});

        tree.computeWeightsOfAncestorsOfLevel(tree.depth);
    }

    if(md && !loadedTree) {
        tree.update();
    }

    if(md || loadedTree) {
        proxy.redraw();
    }
#endif
}

void Viewer::loadExperiment(const char* meshName, const char* treeName, const char* confName) {
    Cr::Utility::Resource expResource{"experiments-data"};

    stanfordImporter.openData(expResource.getRaw(meshName));
    auto md = stanfordImporter.mesh(0);
    mesh.setFromData(*md);
    if(treeName) {
        tree = Tree::deserialize(expResource.getRaw(treeName), mesh);
    } else {
        tree.update();
        tree.root().initializePhasefieldFromParent();
    }
    //tree.update();
    auto raw = expResource.getRaw(confName);
    std::string rawString{raw.begin(), raw.end()};
    std::istringstream is(rawString);
    Cr::Utility::Configuration conf(is);
    loadConfig(conf);
    for(Node node : tree.nodes())
        Debug{} << node;
    Debug{} << "Loaded resources, uploading to gl";
    mesh.uploadVertexBuffer(vertexBuffer);
    mesh.uploadIndexBuffer(indexBuffer);
    glMesh.setCount(mesh.indexCount());
    proxy.redraw();

    setScalingFactors();
    tree.computeWeightsOfAncestorsOfLevel(tree.depth);
}

void Viewer::drawBrushOptions() {
    if(ImGui::TreeNode("Brush")) {
        constexpr int step = 1;
        constexpr double stepDist = 0.01;
        constexpr double min = 0.f, max = 1.f;
        ImGui::SliderScalar("Recursive Phase Filter Factor", ImGuiDataType_Double, &recursiveFilterFactor, &min, &max, "%.3f", ImGuiSliderFlags_Logarithmic);
        ImGui::SliderScalar("Distance Step", ImGuiDataType_Double, &distStep, &min, &max, "%.5f", ImGuiSliderFlags_Logarithmic);
        ImGui::InputScalar("Maximal Distance", ImGuiDataType_Double, &maxDist, &stepDist, nullptr, "%.3f");
        constexpr double lower = -1.f, upper = 1.f;
        ImGui::SliderScalar("Phase", ImGuiDataType_Double, &phase, &lower, &upper, "%.2f", ImGuiSliderFlags_Logarithmic);

        const auto colBrushing = ImVec4(0.56f, 0.83f, 0.26f, 1.0f);
        const auto colNotBrushing = ImVec4(0.85f, 0.85f, 0.85f, 1.0f);
        ImGui::TextColored(brushingModeEnabled ? colBrushing : colNotBrushing, "Press Left Control To Enable Brushing");
        ImGui::SameLine();
        ImVec2 p = ImGui::GetCursorScreenPos();
        ImDrawList* draw_list = ImGui::GetWindowDrawList();
        float height = ImGui::GetFrameHeight();
        float width = height*1.55f;
        ImGui::InvisibleButton("brush modus", ImVec2(width, height));
        draw_list->AddRectFilled(p, ImVec2(p.x + width, p.y + height),
                                 brushingModeEnabled ? ImGui::GetColorU32(colBrushing) : ImGui::GetColorU32(colNotBrushing),
                                 height*0.1f);
        ImGui::TreePop();
    }
}

template<class T>
bool Viewer::drawFunctionals(Array<T>& functionals, size_t& id) {
    bool evaluateProblem = false;
    for(auto& f : functionals) {
        ImGui::PushID(id++);
        ImGui::Separator();

        if(ImGui::Button("Remove")) {
            std::swap(f, functionals.back());
            arrayResize(functionals, functionals.size() - 1);
            ImGui::PopID();
            break;
        }

        ImGui::SameLine();



        if constexpr (std::is_same_v<T, Functional>) {
            if(ImGui::Checkbox("Disable", &f.disable)) {
                proxy.redraw();
            }
            f.drawImGuiOptions(proxy);
        } else {
            if(ImGui::Checkbox("Disable", &f.f.disable)) {
                proxy.redraw();
            }
            f.f.drawImGuiOptions(proxy);
        }

        ImGui::Separator();
        ImGui::PopID();
    }
    return evaluateProblem;
}


void Viewer::drawOptimizationOptions() {
    if(ImGui::TreeNode("Optimization Options")) {

        constexpr static size_t step = 1;

        static auto currentType = FunctionalType::DirichletEnergy;
        ImGui::Text("Functionals:");
        if(ImGui::BeginCombo("##functionals", FunctionalType::to_string(currentType))) {
            for(auto type : FunctionalType::range) {
                bool isSelected = type == currentType;
                if(ImGui::Selectable(FunctionalType::to_string(type), isSelected))
                    currentType = type;
                if(isSelected)
                    ImGui::SetItemDefaultFocus();
            }
            ImGui::EndCombo();
        }

        if(ImGui::Button("Add Functional As Objective"))
            arrayAppend(problem.objectives, {.f = makeFunctional(currentType)});

        if(ImGui::Button("Add Functional As Constraint"))
            arrayAppend(problem.constraints, makeFunctional(currentType));

        bool evaluateProblem = false;

        size_t id = 0;
        evaluateProblem |= drawFunctionals(problem.objectives, id);
        evaluateProblem |= drawFunctionals(problem.constraints, id);

        //static bool drawSegmentation = true;
        //if(ImGui::Checkbox("Segmentation", &drawSegmentation)) {
        //    if(drawSegmentation){
        //        proxy.setTag(segmentationTag);
        //    }
        //}

        //static bool drawPhasefield = true;
        //if(ImGui::Checkbox("Phasefield", &drawPhasefield)) {
        //    if(drawSegmentation){
        //        proxy.setTag(phasefieldTag);
        //    }
        //}

        static bool drawProblemGradient = false;
        if(ImGui::Checkbox("Total Gradient", &drawProblemGradient)) {
            if(drawProblemGradient) {
                proxy.setCallbacks([this](Node) {
                    auto parameters = currentNode.phasefield();
                    Array<double> gradP(parameters.size());
                    double cost = 0;
                    problem.nodeToOptimize = currentNode;
                    totalArea = currentNode.integrateWeight(mesh);
                    problem(parameters, cost, gradP);

                    auto [minimum, maximum] = Math::minmax(gradP);
                    double scale = 1./(maximum - minimum);

                    for(Vertex v : mesh.vertices()) {
                        mesh.scalar(v) = scale*(gradP[v.idx] - minimum);
                    }
                },
                [this] { drawProblemGradient = false; });
            } else {
                proxy.setDefaultCallback();
            }
            proxy.redraw();
        }

        constexpr Double minEps = 0.f, maxEps = 0.1;
        if(ImGui::DragScalar("epsilon", ImGuiDataType_Double, &epsilon, .0001f, &minEps, &maxEps, "%f", 1))
            setScalingFactors();

        constexpr Double sMin = 0, sMax = 5;
        if(ImGui::DragScalar("Kappa", ImGuiDataType_Double, &kappa, .0001f, &sMin, &sMax, "%f", 1))
            setScalingFactors();

        static size_t iterations = 100;
        ImGui::InputScalar("iterations", ImGuiDataType_U64, &options.max_num_iterations, &step, nullptr, "%u");

        auto solverHasSearchDirection = [](Solver::Backend::Value s, Solver::LineSearchDirection::Value direction) {
            if(s == Solver::Backend::IPOPT) {
                switch(direction) {
                    case Solver::LineSearchDirection::LBFGS :
                    case Solver::LineSearchDirection::SR1 :
                        return true;
                    default:
                        return false;
                }
            } else if(s == Solver::Backend::CERES) {
                return direction != Solver::LineSearchDirection::SR1;
            }
            return false;
        };

        if(ImGui::BeginCombo("##solver", Solver::Backend::to_string(options.solver))) {
            for(auto solver : Solver::Backend::range) {
                bool isSelected = (options.solver == solver);
                if(ImGui::Selectable(Solver::Backend::to_string(solver), isSelected)) {
                    options.solver = solver;
                    if(!solverHasSearchDirection(options.solver, options.line_search_direction))
                        options.line_search_direction = Solver::LineSearchDirection::LBFGS;
                }
                if(isSelected)
                    ImGui::SetItemDefaultFocus();
            }
            ImGui::EndCombo();
        }

        if(ImGui::BeginCombo("Descent direction",
                             Solver::LineSearchDirection::to_string(options.line_search_direction))) {
            for(auto dir : Solver::LineSearchDirection::range) {
                if(!solverHasSearchDirection(options.solver, dir)) continue;
                bool isSelected = (options.line_search_direction == dir);
                if(ImGui::Selectable(Solver::LineSearchDirection::to_string(dir), isSelected))
                    options.line_search_direction = dir;
                if(isSelected)
                    ImGui::SetItemDefaultFocus();
            }
            ImGui::EndCombo();
        }

        ImGui::Checkbox("Hierarchical Optimization", &hierarchicalOptimization);

        if(ImGui::Button("Optimize") && !problem.objectives.empty() && !isOptimizing) {
            tree.computeWeightsOfAncestorsOfLevel(currentNode.depth());
            totalArea = currentNode.integrateWeight(mesh);
            isOptimizing = true;

            problem.nodeToOptimize = currentNode;
            pollingSolver.emplace(options, problem, currentNode.phasefield());
        }

        ImGui::SameLine();

        if(ImGui::Button("Stop"))
            isOptimizing = false;

        ImGui::TreePop();

        if(evaluateProblem && !isOptimizing) {
            //Array<Double> gradient(problem.numParameters());
            //double cost = 0;
            //double constraint = 0;
            //problem(tree.phasefieldData, cost, gradient, constraint, nullptr);
        }
    }
}

void Viewer::setCallbacks(UniqueFunction<bool()>&& cb) {
    auto abortCb = [this, cb = std::move(cb)] (Solver::IterationSummary const&) mutable -> Solver::Status::Value {
        proxy.redraw();
        return cb() && isOptimizing ? Solver::Status::CONTINUE : Solver::Status::USER_ABORTED;
    };
    auto& callbacks = options.callbacks;
    arrayAppend(callbacks, InPlaceInit, std::move(abortCb));

    for(auto& o : problem.objectives) arrayResize(o.history, 0);
    t = 0;
    arrayAppend(callbacks, InPlaceInit, [this](Solver::IterationSummary const& summary){
        Node node = problem.nodeToOptimize;
        for(auto& [f, hist, _] : problem.objectives) {
            double cost = 0;
            f(node.phasefield(), node.temporary(), cost, nullptr, nullptr);
            arrayAppend(hist, InPlaceInit, float(t), float(cost));
        }
        ++t;
        return Solver::Status::CONTINUE;
    });
}

void Viewer::runOptimization(UniqueFunction<bool()>&& cb){
    setCallbacks(MOVE(cb));
    if(hierarchicalOptimization) {

        for(Node node : tree.leafs())
            node.splitAndInitialize(&currentNode);
        tree.computeLeafWeights();

        for(Node node : tree.leafs()) {
            totalArea = node.integrateWeight(mesh);
            problem.nodeToOptimize = node;
            Solver::solve(options, problem, node.phasefield(), nullptr);
        }

    } else {
        tree.computeWeightsOfAncestorsOfLevel(currentNode.depth());
        totalArea = currentNode.integrateWeight(mesh);
        problem.nodeToOptimize = currentNode;
        Solver::solve(options, problem, currentNode.phasefield(), nullptr);
    }

    arrayResize(options.callbacks, options.callbacks.size() - 2); /* pop the last two callbacks off */

}

void Viewer::refineLeafNodes(UniqueFunction<bool()>&& cb) {
    setCallbacks(MOVE(cb));

    for(Node leaf : tree.leafs()) {
        leaf.splitAndInitialize(&currentNode);
    }

    tree.computeLeafWeights();

    for(Node leaf : tree.leafs()) {
        totalArea = leaf.integrateWeight(mesh);
        problem.nodeToOptimize = leaf;
        Solver::solve(options, problem, leaf.phasefield(), nullptr);
    }

    arrayResize(options.callbacks, options.callbacks.size() - 2); /* pop the last two callbacks off */
}

Functional Viewer::makeFunctional(FunctionalType::Value type) {
    switch(type) {
        case FunctionalType::AreaRegularizer: {
            AreaRegularizer ar{mesh};
            ar.totalArea = &totalArea;
            Functional f = MOVE(ar);
            f.loss = QuadraticLoss{};
            f.scaling = &areaPenaltyScaling;
            return f;
        }
        case FunctionalType::DirichletEnergy: {
            Functional f = DirichletEnergy{mesh};
            f.scaling = &dirichletScaling;
            return f;
        }
        case FunctionalType::DoubleWellPotential : {
            Functional f = DoubleWellPotential{mesh};
            f.scaling = &doubleWellScaling;
            return f;
        }
        case FunctionalType::ConnectednessConstraint : {
            ConnectednessConstraint cc{mesh};
            cc.epsilon = &epsilon;
            Functional f = MOVE(cc);
            f.scaling = &connectednessScaling;
            return f;
        }
        case FunctionalType::DiffuseYamabe : {
            DiffuseYamabe yamabe{mesh};
            //yamabe.scaling = &yamabeLambdaScaling;
            return yamabe;
        }
    }
    assert(false);
    return Functional{};
}

void Viewer::brush() {
    if(targetDist < maxDist)
        targetDist += distStep;

    Vertex v;
    Double distance;
    while(fastMarchingMethod.step(v, distance)) {
        arrayAppend(distances, InPlaceInit, distance, v);
        if(distance > targetDist)
            break;
    }

    auto phasefield = currentNode.phasefield();
    for(auto [d, v] : distances) {
        auto u = (1.f - recursiveFilterFactor)*phasefield[v.idx] + recursiveFilterFactor*phase;
        phasefield[v.idx] = u;
    }
}

void Viewer::setScalingFactors() {
    dirichletScaling = epsilon/2.;
    connectednessScaling = Math::pow(epsilon, -kappa);
    areaPenaltyScaling = 1.;
    doubleWellScaling = 1./epsilon;
    yamabeLambdaScaling = 1./epsilon;
}

Vector3 Viewer::unproject(Vector2i const& windowPosition, Float depth) {
    /* We have to take window size, not framebuffer size, since the position is
       in window coordinates and the two can be different on HiDPI systems */
    const Vector2i viewSize = windowSize();
    const Vector2i viewPosition{windowPosition.x(), viewSize.y() - windowPosition.y() - 1};
    const Vector3 in{2*Vector2{viewPosition}/Vector2{viewSize} - Vector2{1.0f}, depth*2.0f - 1.0f};

    //get global coordinates
    return (arcBall->transformationMatrix()*projection.inverted()).transformPoint(in);
}

Vertex Viewer::intersectWithPcd(Vector3 const& origin, Vector3 const& dir) {
    BVHAdapter::Intersection intersection;
    size_t idx;
    if(bvh.computeIntersection(origin, dir, intersection)) {
        float u = intersection.u;
        float v = intersection.v;
        size_t vIdx = 0;
        if(u < 0.5 && v > 0.5) {
            vIdx = 1;
        } else {
            vIdx = 2;
        }
        idx =  mesh.triangels()[intersection.idx][vIdx];
        lastIntersectionIdx = idx;
    } else { idx = lastIntersectionIdx; }
    return {idx, &mesh};
}

void Viewer::drawVisualizationOptions() {
    if(ImGui::TreeNode("Visualization Options")) {

        if(ImGui::Button("Redraw")) {
            proxy.redraw();
        }

        if(ImGui::BeginCombo("##visOptions", VisOption::to_string(proxy.option))) {
            for(auto type : VisOption::range) {
                bool isSelected = type == proxy.option;
                if(ImGui::Selectable(VisOption::to_string(type), isSelected)) {
                    proxy.option = type;
                    proxy.setDefaultCallback();
                    proxy.redraw();
                }
                if(isSelected)
                    ImGui::SetItemDefaultFocus();
            }
            ImGui::EndCombo();
        }

        char current[100];
        sprintf(current, "%zu", currentNode.idx);
        if(ImGui::BeginCombo("Current Node", current)) {
            char buffer[100];
            for(Node node : tree.nodes()) {
                bool isSelected = node == currentNode;
                sprintf(buffer, "%zu", node.idx);
                if(ImGui::Selectable(buffer, isSelected)) {
                    currentNode = node;
                    if(proxy.option != VisOption::Segmentation)
                        proxy.redraw();
                }
                if(isSelected)
                    ImGui::SetItemDefaultFocus();
            }
            ImGui::EndCombo();
        }

        if(ImGui::Button("Optimize Colors")) {
            Array<size_t> starts, neighbors;
            tree.computeAdjacencyGraph(neighbors, starts, proxy.level);
            Debug{} << getColors();
            optimizeColors(neighbors, starts);
            Debug{} << getColors();
            if(proxy.option == VisOption::Segmentation)
                proxy.redraw();

        }

        if(ImGui::Checkbox("Swap Colors", &swapColors)) {
            if(!swapColors) {
                swapIndex = 0;
            }
        }

        if(ImGui::InputInt("Segmentation level", &proxy.level)) {
            proxy.level = Math::clamp<int>(proxy.level, 0, tree.depth);
            if(proxy.option == VisOption::Segmentation)
                proxy.redraw();
        }

        if(ImGui::Button("Split Leaf nodes")) {
            for(Node node : tree.leafs()) {
                node.splitAndInitialize(&currentNode);
            }
            proxy.redraw();
        }

        if(ImGui::Button("Initialize Current Node")) {
            currentNode.initializePhasefieldFromParent();
            proxy.redraw();
        }

        if(ImGui::Button("Split Current Node")) {
            /* the node handle (potentially) gets invalidated after we add the other child */
            Node rightChild = currentNode.addRightChild(&currentNode);
            rightChild.initializePhasefieldFromParent();

            Node leftChild = currentNode.addLeftChild(&currentNode);
            leftChild.initializePhasefieldFromParent();
            proxy.redraw();
        }

        if(ImGui::Button("Reset Phasefield Tree")) {
            tree.reset();
            currentNode = tree.root();
            proxy.redraw();
        }

        if(ImGui::Button("Reset Current Phase")) {
            for(double& x : currentNode.phasefield()) x = 0;
            proxy.redraw();
        }

        if(ImGui::Button("Print Information")) {

            DirichletEnergy de(mesh);
            DoubleWellPotential dw(mesh);

            SmootherStep chi;
            Array<Pair<double, double>> patchInformation;
            tree.computeLeafWeights();
            for(Node leaf : tree.leafs()) {
                auto phasefield = leaf.phasefield();
                auto weights = leaf.temporary();
                double posArea = 0;
                double negArea = 0;
                for(Face f : mesh.faces()) {
                    double x = 0;
                    double y = 0;
                    for(Vertex v : f.vertices()) {
                        x += chi.eval(phasefield[v])*weights[v];
                        y += chi.eval(-phasefield[v])*weights[v];
                    }
                    posArea += f.area()*x/3;
                    negArea += f.area()*y/3;
                }

                double deResult = 0, dwResult = 0;
                de.operator()<double>(phasefield, weights, deResult, nullptr, nullptr);
                dw.operator()<double>(phasefield, weights, dwResult, nullptr, nullptr);

                double boundary = doubleWellScaling*dwResult+dirichletScaling*deResult;
                arrayAppend(patchInformation, {{posArea, boundary}, {negArea, boundary}});
            }

            CORRADE_ASSERT(patchInformation.size() == tree.numLeafs*2, "weird number of patches",);
            double targetArea = 0;
            for(auto [x,_] : patchInformation) targetArea += x;
            targetArea /= double(patchInformation.size());
            printf("Target Area %f\n", targetArea);
            double totalError = 0;

            auto& colors = getColors(tree.numLeafs*2);
            CORRADE_ASSERT(colors.size() == patchInformation.size(), "Patch areas not same length as colors",);

            for(size_t i = 0; i < colors.size(); ++i) {
                auto [patchArea, patchBoundary] = patchInformation[i];
                double error = Math::abs(targetArea - patchArea);
                totalError += error;
                printf("Patch Area Error %f (area = %f), Boundary Length %f\n", error, patchArea, patchBoundary);
                Debug{} << "Patch Color" << colors[i];
            }
            printf("Total error %f\n", totalError);
        }

        static bool drawCurvature = false;
        if(ImGui::Checkbox("Gaussian Curvature", &drawCurvature)) {
            if(drawCurvature) {
                proxy.setCallbacks(
                        [this](Node node) {
                            mesh.requireGaussianCurvature();
                            proxy.drawValuesNormalized(mesh.gaussianCurvature);
                        },
                        [this]{ drawCurvature = false; });
            } else proxy.setDefaultCallback();
            proxy.redraw();
        }

        if(ImGui::Checkbox("Custom Mapping", &proxy.customMapping)) {
            proxy.redraw();
        }

        static const double minScaling = 0.00001, maxScaling = 10;
        static const double minOffset = -1, maxOffset = 1;

        if(ImGui::DragScalar("Scaling", ImGuiDataType_Double, &proxy.scale, .001f, &minScaling, &maxScaling, "%f"))
            proxy.redraw();

        ImGui::SameLine();

        if(ImGui::DragScalar("Offset", ImGuiDataType_Double, &proxy.offset, .001f, &minOffset, &maxOffset, "%f"))
            proxy.redraw();


        if(ImGui::BeginCombo("Color Map", colorMapData[colorMapIndex].name)) {
            for(size_t i = 0; i < colorMapData.size(); ++i) {
                bool isSelected = i == colorMapIndex;
                if(ImGui::Selectable(colorMapData[i].name, isSelected)) {
                    colorMapIndex = i;
                    if(proxy.shaderConfig == VisualizationProxy::ShaderConfig::ColorMaps)
                        redraw();
                }
                if(isSelected)
                    ImGui::SetItemDefaultFocus();
            }
            ImGui::EndCombo();
        }

        if(ImGui::Checkbox("Draw Wireframe", &drawWireFrame)) {

        }


        if(ImGui::Button("Bake To Vertex Color")) {
            for(Vertex v : mesh.vertices()) {
                float coord = mesh.scalar(v);
                auto idx = size_t(coord*255);
                mesh.color(v) = Math::unpack<Color3>(colorMapData[colorMapIndex].colors[idx]);
            }
        }

        if(ImGui::Button("Show Error Plot")) {
            showPlot = true;
        }

        ImGui::TreePop();
    }
}

void Viewer::drawIO() {
    if(ImGui::TreeNode("IO")) {


        static size_t curExp = 0;

        if(ImGui::BeginCombo("Experiment", experiments[curExp].name)) {
            for(size_t i = 0; i < IM_ARRAYSIZE(experiments); ++i) {
                bool isSelected = i == curExp;
                if(ImGui::Selectable(experiments[i].name, isSelected)) {
                    curExp = i;
                }
                if(isSelected)
                    ImGui::SetItemDefaultFocus();
            }
            ImGui::EndCombo();
        }

        if(ImGui::Button("Load Experiment")) {
            auto& exp = experiments[curExp];
            loadExperiment(exp.meshName, exp.treeName, exp.confName);
        }


        static char path[50] = "/home/janos/";
        static char postfix[25] = "torus_hierarchical_32";

        ImGui::BeginGroup();
        ImGui::Text("Import/Export to");
        ImGui::InputText("##Export to", path, sizeof(path));
        ImGui::EndGroup();

        ImGui::SameLine();

        ImGui::BeginGroup();
        ImGui::Text("Mesh Name");
        ImGui::InputText("##Mesh Name", postfix, sizeof(postfix));
        ImGui::EndGroup();

        ImGui::Dummy(ImVec2(0.0f, 20.0f));

        if(ImGui::Button("Export Scene")) {
            if(!Cr::Utility::Directory::exists(path))
                Cr::Utility::Directory::mkpath(path);

            Array<char> data;
            tree.serialize(data);

            if(strlen(postfix)) {
                std::string postfixStr{postfix};
                std::string treePath = Cr::Utility::Directory::join(path, postfixStr + ".bin");
                FILE* fp = fopen(treePath.c_str(), "w");
                fwrite(data.data(), sizeof(char), data.size(), fp);
                fclose(fp);

                std::string meshPath = Cr::Utility::Directory::join(path, postfixStr + ".ply");
                saveMesh(meshPath.c_str());
                std::string confPath = Cr::Utility::Directory::join(path, postfixStr + ".conf");
                saveCurrentConfig(confPath.c_str());
            }
        }

        ImGui::SameLine();

        if(ImGui::Button("Load Scene")) {
            loadScene(path, postfix);
        }

        if(ImGui::Button("Export to VTK")) {
            if(!Cr::Utility::Directory::exists(path))
                Cr::Utility::Directory::mkpath(path);

            std::string meshPath = Cr::Utility::Directory::join(path,"dump.vts");
            dumpMesh(meshPath.c_str());
        }

        ImGui::Checkbox("Animate", &animate);

        static char recordingPath[100] = "/home/janos/test.mp4";
        ImGui::InputText("Video Path", recordingPath, 100);

        const auto green = ImVec4(0.56f, 0.83f, 0.26f, 1.0f);
        const auto red = ImVec4(0.83f, 0.56f, 0.26f, 1.0f);

        ImGui::PushStyleColor(ImGuiCol{}, recording ? red : green);
        if(ImGui::Button("Start Recording")) {
            if(!recording) {
#ifdef PHASEFIELD_WITH_IO
                videoSaver.startRecording(recordingPath, framebufferSize());
                recording = true;
#endif
            }
        }
        ImGui::PopStyleColor();
        ImGui::SameLine();
        ImGui::PushStyleColor(ImGuiCol{}, recording ? green : red);
        if(ImGui::Button("Stop Recording")) {
            if(recording) {
#ifdef PHASEFIELD_WITH_IO
                videoSaver.endRecording();
                recording = false;
#endif
            }
        }
        ImGui::PopStyleColor();

        ImGui::TreePop();
    }
}

void Viewer::viewportEvent(ViewportEvent& event) {
    Mg::GL::defaultFramebuffer.setViewport({{}, event.framebufferSize()});
    arcBall->reshape(event.windowSize());
    projection = Matrix4::perspectiveProjection(fov, Vector2{framebufferSize()}.aspectRatio(), 0.01f, 100.0f);

    imgui.relayout(Vector2{event.windowSize()}/event.dpiScaling(),
                   event.windowSize(), event.framebufferSize());

}

void Viewer::keyPressEvent(KeyEvent& event) {
    if(imgui.handleKeyPressEvent(event)) {
        event.setAccepted();
        return;
    }

    switch(event.key()) {
        case KeyEvent::Key::R:
            arcBall->reset();
            break;
        case KeyEvent::Key::LeftCtrl :
            brushingModeEnabled = true;
            break;
        default:
            return;
    }

    event.setAccepted();
    redraw(); /* camera or mesh has changed, redraw! */
}

void Viewer::keyReleaseEvent(KeyEvent& event) {
    if(imgui.handleKeyReleaseEvent(event)) {
        event.setAccepted();
        return;
    }

    if(event.key() == Viewer::KeyEvent::Key::LeftCtrl) {
        brushingModeEnabled = false;
        brushing = false;
        event.setAccepted();
        return;
    }
}

void Viewer::textInputEvent(TextInputEvent& event) {
    if(imgui.handleTextInputEvent(event)) {
        event.setAccepted();
        return;
    }
}

void Viewer::startBrushing(Vector3 const& origin, Vector3 const& dir) {
    fastMarchingMethod.update();

    arrayResize(distances, 0);
    Vertex v = intersectWithPcd(origin, dir);
    fastMarchingMethod.setSource(v);
    targetDist = 0.;
    brushing = true;
}

void Viewer::mousePressEvent(MouseEvent& event) {
    if(imgui.handleMousePressEvent(event)) {
        event.setAccepted();
        return;
    }

    if(brushingModeEnabled) {
        Vector3 o = unproject(event.position(), 0);
        Vector3 d = unproject(event.position(), 1) - o;
        startBrushing(o, d);
        event.setAccepted();
        return;
    }

    if(swapColors) {
        Vector3 o = unproject(event.position(), 0);
        Vector3 d = unproject(event.position(), 1) - o;
        Vertex v = intersectWithPcd(o, d);
        Color4 target = mesh.color(v);
        auto& colors = getColors();
        float min = std::numeric_limits<float>::max();
        size_t idx = Invalid;
        for(size_t i = 0; i < colors.size(); ++i) {
            float dist = (colors[i] - target).dot();
            if(dist < min) {
                min = dist;
                idx = i;
            }
        }
        colorIndexToSwap[swapIndex] = idx;
        if(swapIndex) {
            std::swap(colors[colorIndexToSwap[0]], colors[colorIndexToSwap[1]]);
            proxy.redraw();
        }
        swapIndex ^= 1;
        event.setAccepted();
        return;
    }

    if(event.button() == MouseEvent::Button::Right) {
        Vector3 o = unproject(event.position(), 0);
        Vector3 d = unproject(event.position(), 1) - o;

        Vertex v = intersectWithPcd(o, d);
        Debug{} << "Phasefield value at mouse location" << currentNode.phasefield()[v];
        event.setAccepted();
        return;
    }

    trackingMouse = true;

    arcBall->initTransformation(event.position());

    event.setAccepted();
    redraw(); /* camera has changed, redraw! */
}

void Viewer::mouseReleaseEvent(MouseEvent& event) {
    if(imgui.handleMouseReleaseEvent(event)) {
        event.setAccepted();
        return;
    }

    if(brushingModeEnabled) {
        brushing = false;
        event.setAccepted();
        return;
    }

    /* Disable mouse capture again */
    /** @todo replace once https://github.com/mosra/magnum/pull/419 is in */
    if(trackingMouse) {

        trackingMouse = false;
        event.setAccepted();
    }
}

void Viewer::mouseMoveEvent(MouseMoveEvent& event) {
    if(imgui.handleMouseMoveEvent(event)) {
        event.setAccepted();
        return;
    }

    if(brushingModeEnabled && brushing){
        Vector3 o = unproject(event.position(), 0);
        Vector3 d = unproject(event.position(), 1) - o;
        startBrushing(o, d);
        event.setAccepted();
        return;
    }

    if(trackingMouse) {
        if(event.modifiers() & MouseMoveEvent::Modifier::Shift)
            arcBall->translate(event.position());
        else arcBall->rotate(event.position());

        event.setAccepted();
        redraw(); /* camera has changed, redraw! */
    }
}

void Viewer::mouseScrollEvent(MouseScrollEvent& event) {
    if(imgui.handleMouseScrollEvent(event)) {
        /* Prevent scrolling the page */
        event.setAccepted();
        return;
    }

    const Float delta = event.offset().y();
    if(Math::abs(delta) < 1.0e-2f) return;

    arcBall->zoom(delta);

    event.setAccepted();
    redraw(); /* camera has changed, redraw! */
}

void Viewer::drawEvent() {
    Mg::GL::defaultFramebuffer.clear(Mg::GL::FramebufferClear::Color | Mg::GL::FramebufferClear::Depth);
    imgui.newFrame();

    /* Enable text input, if needed */
    if(ImGui::GetIO().WantTextInput && !isTextInputActive())
        startTextInput();
    else if(!ImGui::GetIO().WantTextInput && isTextInputActive())
        stopTextInput();

    if(brushing) {
        brush();
        proxy.redraw();
    }

    if(animate) {
        Vector2i m = windowSize()/2;
        arcBall->initTransformation(m);
        arcBall->rotate(m + Vector2i{5, 0});
    }

#ifdef MAGNUM_TARGET_WEBGL
    if(isOptimizing) {
        if(pollingSolver->runOneIteration() < 0)
            isOptimizing = false;
        proxy.redraw();
    }
#endif

    proxy.upload(); /* synchronize with gpu */

    /* draw scene */
    bool camChanged = arcBall->updateTransformation();
    Matrix4 viewTf = arcBall->viewMatrix();

    if(proxy.shaderConfig == VisualizationProxy::ShaderConfig::ColorMaps) {
        phongColorMap.setProjectionMatrix(projection)
                     .setTransformationMatrix(viewTf)
                     .setNormalMatrix(viewTf.normalMatrix())
                     .bindDiffuseTexture(colorMapData[colorMapIndex].texture)
                     .draw(glMesh);
    } else if (proxy.shaderConfig == VisualizationProxy::ShaderConfig::VertexColors) {
        phongVertexColors.setProjectionMatrix(projection)
                         .setTransformationMatrix(viewTf)
                         .setNormalMatrix(viewTf.normalMatrix())
                         .setLightPositions({{10,10,10, 0}, {-10, -10, 10, 0}})
                         .setLightColors({Color3{0.5}, Color3{0.5}})
                         .setSpecularColor(Color4{0.1})
                         .draw(glMesh);
    }

    if(drawWireFrame) {
        GL::Renderer::setDepthFunction(GL::Renderer::DepthFunction::LessOrEqual);
        GL::Renderer::enable(GL::Renderer::Feature::Blending);
        //GL::Renderer::setBlendFunction(
        //        GL::Renderer::BlendFunction::One, /* or SourceAlpha for non-premultiplied */
        //        GL::Renderer::BlendFunction::OneMinusSourceAlpha);

        StridedArrayView1D<const UnsignedInt> indices = mesh.indices();
        StridedArrayView1D<const Vector3> indexedPositions = mesh.positions();

        /* De-indexing the position array */
        GL::Buffer vertices{Mg::MeshTools::duplicate(indices, indexedPositions)};

        GL::Mesh wireframeMesh;
        wireframeMesh.addVertexBuffer(std::move(vertices), 0, Mg::Shaders::MeshVisualizer3D::Position{})
                     .setCount(indices.size());

        meshVis.setColor(Color4{0,0,0,0})
               .setWireframeColor(0xdcdcdc_rgbf)
               .setTransformationMatrix(viewTf)
               .setProjectionMatrix(projection)
               .draw(wireframeMesh);

        GL::Renderer::disable(GL::Renderer::Feature::Blending);
        GL::Renderer::setDepthFunction(GL::Renderer::DepthFunction::Less);
    }


    if(recording) {
#ifdef PHASEFIELD_WITH_IO
        Mg::Image2D image = GL::defaultFramebuffer.read({{},framebufferSize()}, {GL::PixelFormat::RGBA, Mg::GL::PixelType::UnsignedByte});
        videoSaver.appendFrame(std::move(image));
#endif
    }


    /* draw ImGui stuff */

    ImGui::Begin("Phasefield Options");
    ImGui::PushItemWidth(150);
    drawBrushOptions();
    drawOptimizationOptions();
    drawVisualizationOptions();
    drawIO();
    ImGui::PopItemWidth();
    ImGui::End();

    if(showPlot)
        drawErrorPlot();

    imgui.updateApplicationCursor(*this);

    /* Render ImGui window */
    {
        Mg::GL::Renderer::enable(Mg::GL::Renderer::Feature::Blending);
        Mg::GL::Renderer::disable(Mg::GL::Renderer::Feature::FaceCulling);
        Mg::GL::Renderer::disable(Mg::GL::Renderer::Feature::DepthTest);
        Mg::GL::Renderer::enable(Mg::GL::Renderer::Feature::ScissorTest);

        imgui.drawFrame();

        Mg::GL::Renderer::disable(Mg::GL::Renderer::Feature::ScissorTest);
        Mg::GL::Renderer::enable(Mg::GL::Renderer::Feature::DepthTest);
        Mg::GL::Renderer::enable(Mg::GL::Renderer::Feature::FaceCulling);
        Mg::GL::Renderer::disable(Mg::GL::Renderer::Feature::Blending);
    }


    swapBuffers();
    redraw();
}

Viewer::~Viewer() {
    ImPlot::DestroyContext();
}

void Viewer::drawErrorPlot() {
    ImGui::Begin("Optimization Error", &showPlot, ImGuiWindowFlags_MenuBar);

    size_t objectiveCount = problem.objectives.size();

    if (ImGui::Button("Clear", ImVec2(100,0))) {
        t = 0;
        for (auto& [f, history, show] : problem.objectives)
            arrayResize(history, 0);
    }

    ImGui::SameLine();

    if(ImGui::Button("Save Tex plot", ImVec2(100, 0))) {
        std::string plot = "\\begin{tikzpicture}\n"
                           "\\begin{axis}[\n"
                           "height=9cm,\n"
                           "width=9cm,\n"
                           "grid=major,\n"
                           "]\n";

        size_t step = Math::max(t/50, 1ul);
        for(auto& [f,hist,show] : problem.objectives) {
            plot += "\n\\addplot coordinates {\n";
            for(size_t j = 0; j < hist.size(); j += step) {
                plot += Cr::Utility::formatString("({},{})\n", j, hist[j].y());
            }
            plot += "};\n";
            plot += Cr::Utility::formatString("\\addlegendentry{{ {} }}",
                                              FunctionalType::to_string(f.functionalType));
        }
        plot += "\n\\end{axis}\n"
                "\\end{tikzpicture}";

        FILE* fp = fopen("/tmp/plot.tex", "w");
        if(fp != nullptr) {
            fputs(plot.data(), fp);
            fclose(fp);
        }
    }


    ImPlot::SetNextPlotLimitsX(0, t, paused ? ImGuiCond_Once : ImGuiCond_Always);
    if (ImPlot::BeginPlot("##DND", nullptr, nullptr, ImVec2(-1,0), ImPlotFlags_YAxis2 | ImPlotFlags_YAxis3, ImPlotAxisFlags_NoTickLabels)) {
        for (auto& [f, history, show] : problem.objectives) {
            const char* label = FunctionalType::to_string(f.functionalType);
            if (show && history.size() > 0) {
                ImPlot::SetPlotYAxis(0);
                ImPlot::PlotLine(label, &history[0].x(), &history[0].y(), history.size(), 0, 2 * sizeof(float));
            }
        }

        ImPlot::EndPlot();
    }
    ImGui::End();
}

#ifdef MAGNUM_TARGET_WEBGL
Int Viewer::touchStartEvent(EmscriptenTouchEvent const* event) {
    if (event->numTouches == 2) {
        isPinching = true;
        EmscriptenTouchPoint p1 = event->touches[0];
        EmscriptenTouchPoint p2 = event->touches[1];
        Vector2d p1Client{double(p1.clientX), double(p1.clientY)};
        Vector2d p2Client{double(p2.clientX), double(p2.clientY)};
        pinchLength = (p1Client - p2Client).length();
        return 1;
    } else if(event->numTouches == 1) {
        EmscriptenTouchPoint ep = event->touches[0];
        Vector2i p{Int(ep.targetX), Int(ep.targetY)};

        ImGuiIO& io = ImGui::GetIO();
        io.MousePos = ImVec2(dpiScaling()*Vector2(p));
        io.MouseDown[0] = true;
        if(io.WantCaptureMouse) {
            trackingForImGui = true;
            redraw();
            return true;
        }

        arcBall->initTransformation(p);
        trackingFinger = true;
        return 1;
    } else if(event->numTouches == 3) {
        EmscriptenTouchPoint ep = event->touches[0];
        arcBall->initTransformation({Int(ep.targetX), Int(ep.targetY)});
        trackingFingers = true;
        return 1;
    }
    return 0;
}

Int Viewer::touchMoveEvent(EmscriptenTouchEvent const* event) {
    ImGuiIO& io = ImGui::GetIO();

    if(trackingForImGui && io.WantCaptureMouse) {
        EmscriptenTouchPoint ep = event->touches[0];
        Vector2 p{Float(ep.targetX), Float(ep.targetY)};
        io.MousePos = ImVec2(Vector2(p));
        return 1;
    }

    if (isPinching) {
        EmscriptenTouchPoint p1 = event->touches[0];
        EmscriptenTouchPoint p2 = event->touches[1];
        Vector2d p1Client{double(p1.targetX), double(p1.targetY)};
        Vector2d p2Client{double(p2.targetY), double(p2.targetY)};
        double length = (p1Client - p2Client).length();
        double delta = 5.*(length - pinchLength)/Float(windowSize().max());
        arcBall->zoom(delta);
        pinchLength = length;
        return 1;
    } else if (trackingFinger) {
        EmscriptenTouchPoint ep = event->touches[0];
        arcBall->rotate({Int(ep.targetX), Int(ep.targetY)});
        return 1;
    } else if(trackingFingers) {
        EmscriptenTouchPoint ep = event->touches[0];
        arcBall->translate({Int(ep.targetX), Int(ep.targetY)});
        return 1;
    }
    return 0;
}

Int Viewer::touchEndEvent(EmscriptenTouchEvent const* event) {
    if(trackingForImGui) {
        trackingForImGui = false;
        ImGuiIO& io = ImGui::GetIO();
        io.MouseDown[0] = false;
        return 1;
    }
    if(isPinching || trackingFinger || trackingFingers) {
        isPinching = false;
        trackingFinger = false;
        trackingFingers = false;
        return 1;
    }
    return 0;
}


Int Viewer::touchCancelEvent(EmscriptenTouchEvent const* event) {
    if(trackingForImGui) {
        trackingForImGui = false;
        ImGuiIO& io = ImGui::GetIO();
        io.MouseDown[0] = false;
        return 1;
    }
    if(isPinching || trackingFinger || trackingFingers) {
        isPinching = false;
        trackingFinger = false;
        trackingFingers = false;
        return 1;
    }
    return 0;
}

void Viewer::drawMeshEdit() {
    if(ImGui::TreeNode("Mesh Editor")) {




        ImGui::TreePop();
    }
}
#endif


}
