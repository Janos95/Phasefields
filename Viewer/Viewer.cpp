//
// Created by janos on 08.11.19.
//

#include "Viewer.h"
#include "ImGuiWidgets.h"
#include "Enums.h"
#include "Tag.h"
#include "ModicaMortola.h"
#include "ConnectednessConstraint.h"
#include "WeakYamabe.h"

#include <ScopedTimer/ScopedTimer.h>

#include <Corrade/Utility/Algorithms.h>
#include <Corrade/Utility/Configuration.h>
#include <Corrade/Utility/Resource.h>
#include <Corrade/Utility/Directory.h>
#include <Corrade/PluginManager/Manager.h>
#include <Corrade/PluginManager/PluginMetadata.h>
#include <Corrade/Containers/GrowableArray.h>
#include <Corrade/Utility/ConfigurationGroup.h>

#include <Magnum/GL/DefaultFramebuffer.h>
#include <Magnum/GL/Renderer.h>
#include <Magnum/GL/PixelFormat.h>
#include <Magnum/GL/TextureFormat.h>
#include <Magnum/GL/Buffer.h>
#include <Magnum/GL/Mesh.h>
#include <Magnum/Trade/AbstractSceneConverter.h>

#include <Magnum/Math/Matrix4.h>
#include <Magnum/Math/FunctionsBatch.h>

#include <Magnum/Shaders/Flat.h>
#include <Magnum/Shaders/MeshVisualizer.h>
#include <Magnum/Shaders/Phong.h>
#include <Magnum/Shaders/VertexColor.h>

#include <Magnum/Trade/MeshData.h>
#include <Magnum/Trade/AbstractImporter.h>
#include <Magnum/Image.h>
#include <Magnum/Primitives/Axis.h>
#include <Magnum/ImageView.h>
#include <Magnum/PixelFormat.h>
#include <Magnum/DebugTools/ColorMap.h>
#include <Magnum/DebugTools/ObjectRenderer.h>
#include <Magnum/ImGuiIntegration/Context.hpp>

#include <Magnum/MeshTools/Compile.h>
#include <Magnum/MeshTools/RemoveDuplicates.h>

#include <Magnum/Primitives/Capsule.h>
#include <Magnum/Primitives/Cylinder.h>
#include <Magnum/Primitives/Grid.h>
#include <MagnumPlugins/PrimitiveImporter/PrimitiveImporter.h>

namespace Phasefield {

using namespace Mg::Math::Literals;

namespace {

bool loadMesh(char const* path, Mg::Trade::MeshData& mesh) {
    Mg::PluginManager::Manager<Mg::Trade::AbstractImporter> manager;
    auto importer = manager.loadAndInstantiate("AnySceneImporter");
    auto name = importer->metadata()->name();
    Debug{} << "Trying to load mesh using " << name.c_str();
    if(!importer) std::exit(1);

    Debug{} << "Opening file" << path;

    if(!importer->openFile(path)) {
        puts("could not open file");
        std::exit(4);
    }

    Debug{} << "Imported " << importer->meshCount() << " meshes";

    if(importer->meshCount() && importer->mesh(0)) {
        mesh = *importer->mesh(0);
        return true;
    } else {
        Debug{} << "Could not load mesh";
        return false;
    }
}



bool saveMesh(char const* path, Mesh const& mesh) {
    Mg::PluginManager::Manager<Mg::Trade::AbstractSceneConverter> manager;
    auto exporter = manager.loadAndInstantiate("AnySceneConverter");

    Mg::Trade::MeshData md = mesh.meshDataView();
    if(!exporter->convertToFile(path, md)) {
        Mg::Error{} << "Cannot save file to " << path;
        return false;
    }
    return true;
}

Array<std::pair<char const*, Mg::GL::Texture2D>> makeColorMapTextures() {

    Array<std::pair<char const*, Mg::GL::Texture2D>> textures;
    using L = std::initializer_list<std::pair<char const*, StaticArrayView<256, const Vector3ub>>>;
    for(auto&& [name, colorMap] : L{
            {"Turbo",   Magnum::DebugTools::ColorMap::turbo()},
            {"Magma",   Magnum::DebugTools::ColorMap::magma()},
            {"Plasma",  Magnum::DebugTools::ColorMap::plasma()},
            {"Inferno", Magnum::DebugTools::ColorMap::inferno()},
            {"Viridis", Magnum::DebugTools::ColorMap::viridis()}
    }) {
        const Magnum::Vector2i size{Magnum::Int(colorMap.size()), 1};
        Mg::GL::Texture2D texture;
        texture.setMinificationFilter(Magnum::SamplerFilter::Linear)
               .setMagnificationFilter(Magnum::SamplerFilter::Linear)
               .setWrapping(Magnum::SamplerWrapping::ClampToEdge) // or Repeat
               .setStorage(1, Mg::GL::TextureFormat::RGB8, size) // or SRGB8
               .setSubImage(0, {}, Mg::ImageView2D{Magnum::PixelFormat::RGB8Srgb, size, colorMap});
        arrayAppend(textures, InPlaceInit, name, std::move(texture));
    }
    return textures;
}

}

//Solver::Status::Value OptimizationCallback::operator()(Solver::IterationSummary const&) {
//    if(!optimize) return Solver::Status::USER_ABORTED;
//    return Solver::Status::CONTINUE;
//}


Viewer::Viewer(int argc, char** argv) :
        Mg::Platform::Application{{argc, argv}, Mg::NoCreate},
        fastMarchingMethod(mesh),
        tree{mesh},
        kdtree{mesh},
        proxy(*this),
        dirichletScaling(1.), doubleWellScaling(1.), connectednessScaling(1.),
        problem(tree)
        //problem(tree),
        //proxy(*this),
        //segmentationTag(getTag()),
        //phasefieldTag(getTag())
{
    {
        //arrayAppend(options.callbacks, InPlaceInit, optimizationCallback);
        currentNode = Node{0, &tree};
    }

    /* Setup window */
    {
        ScopedTimer t{"Creating OpenGL Context", true};
        const Vector2 dpiScaling = this->dpiScaling({});
        Configuration conf;
        conf.setTitle("Viewer")
            .setSize({1200, 1200}, dpiScaling)
            .setWindowFlags(Configuration::WindowFlag::Resizable);
        GLConfiguration glConf;
        glConf.setSampleCount(dpiScaling.max() < 2.0f ? 8 : 2);
        if(!tryCreate(conf, glConf)) {
            create(conf, glConf.setSampleCount(0));
        }
    }

    /* setup shaders, color map textures and mesh*/
    {
        ScopedTimer t{"Setting up opengl stuff (textures etc)", true};
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
        vertexBuffer = Mg::GL::Buffer{};
        indexBuffer = Mg::GL::Buffer{};

        mesh.uploadVertexBuffer(vertexBuffer);
        mesh.uploadIndexBuffer(indexBuffer);

        glMesh.setPrimitive(Mg::MeshPrimitive::Triangles)
              .setCount(mesh.indexCount())
              .setIndexBuffer(indexBuffer, 0, Mg::MeshIndexType::UnsignedInt)
              .addVertexBuffer(vertexBuffer, 0, Phong::Position{}, Phong::Normal{}, Phong::TextureCoordinates{}, Phong::Color4{});

        phong = Phong{Phong::Flag::VertexColor, 2};

        colorMapTextures = makeColorMapTextures();

    }

    /* setup the problem */
    {
        arrayAppend(problem.objectives, makeFunctional(FunctionalType::DoubleWellPotential));
        arrayAppend(problem.objectives, makeFunctional(FunctionalType::DirichletEnergy));
        arrayAppend(problem.objectives, makeFunctional(FunctionalType::AreaRegularizer));
        arrayAppend(problem.objectives, makeFunctional(FunctionalType::ConnectednessConstraint));
    }

    /* try to load tree from disk, otherwise fallback to empty tree with three nodes */
    {

        Mg::PluginManager::Manager<Mg::Trade::AbstractImporter> manager;
        primitiveImporter = new Mg::Trade::PrimitiveImporter{manager, "PrimitiveImporter"};
        primitiveImporter->openData("deadbeef");

        //primitiveImporter = manager.instantiate("PrimitiveImporter");

        ScopedTimer t{"Loading scene from disk", true};
        loadScene("/home/janos/meshes/testing", "spiky_cube");

    }

    /* Setup ImGui, load a better font */
    {
        ImGui::CreateContext();
        ImGui::StyleColorsDark();

        //ImFontConfig fontConfig;
        //fontConfig.FontDataOwnedByAtlas = false;
        //const Vector2 size = Vector2{windowSize()}/dpiScaling();
        //Cr::Utility::Resource rs{"data"};
        //ArrayView<const char> font = rs.getRaw("SourceSansPro-Regular.ttf");
        //ImGui::GetIO().Fonts->AddFontFromMemoryTTF(
        //        const_cast<char*>(font.data()), Int(font.size()), 16.0f*framebufferSize().x()/size.x(), &fontConfig);

        imgui = Mg::ImGuiIntegration::Context(*ImGui::GetCurrentContext(),
                                                  Vector2{windowSize()}/dpiScaling(), windowSize(), framebufferSize());

        /* Setup proper blending to be used by ImGui */
        GL::Renderer::setBlendEquation(
                GL::Renderer::BlendEquation::Add, GL::Renderer::BlendEquation::Add);
        GL::Renderer::setBlendFunction(
                GL::Renderer::BlendFunction::SourceAlpha,
                GL::Renderer::BlendFunction::OneMinusSourceAlpha);
    }

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

    /* Start the timer, loop at 60 Hz max */
    setSwapInterval(1);
    setMinimalLoopPeriod(16);
}

void Viewer::loadScene(const char* path, const char* postfix) {
    bool loadedTree = false;
    Optional<Mg::Trade::MeshData> md;

    if(Cr::Utility::Directory::exists(path)) {

        std::string meshPath = Cr::Utility::Directory::join(path, std::string{postfix} + ".ply");
        std::string treePath = Cr::Utility::Directory::join(path, std::string{postfix} + ".bin");


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

            Debug{} << "Opening file" << meshPath;

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

    } else { /* try a primitive */
        Cr::Utility::Configuration conf{"/home/janos/Phasefield/primitives.conf"};
        auto* group = conf.group(postfix);
        auto* importerGroup = primitiveImporter->configuration().group(postfix);

        //primitiveImporter->configuration() = *conf.group("configuration");
        md = primitiveImporter->mesh(postfix);
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
}

//template<class T>
//Array<Vector3d> positionsAsArray(Mg::Trade::MeshData const& meshData) {
//    auto view = meshData.attribute<Vector3>(Mg::Trade::MeshAttribute::Position);
//    Array<Vector3d> positions(NoInit, view.size());
//    for(int i = 0; i < view.size(); ++i) {
//        positions[i] = Vector3d{view[i]};
//    }
//    return positions;
//}

//void Viewer::drawSubdivisionOptions() {
//    if(ImGui::TreeNode("Subdivisions")) {
//        ImGui::Text("Currently we have %d vertices and %d faces", (int) vertices.size(), (int) indices.size()/3);
//        constexpr int step = 1;
//        auto old = numSubdivisions;
//        if(ImGui::InputScalar("Number of Sibdivision (wrt. orignal mesh)", ImGuiDataType_U32, &numSubdivisions, &step,
//                              nullptr, "%d")) {
//            if(numSubdivisions > old) {
//                tree.subdivide(indices, vertices);
//            } else {
//                indices = original.indicesAsArray();
//                vertices = positionsAsArray<Double>(original);
//                tree.resize(vertices.size());
//                auto copy = numSubdivisions;
//                while(copy--)
//                    tree.subdivide(indices, vertices);
//            }
//
//            updateInternalDataStructures();
//            upload(mesh, vertexBuffer, indexBuffer, meshData);
//        }
//
//        ImGui::TreePop();
//    }
//}

//void normalizeMesh(Array<Vector3d>& vertices, Array<UnsignedInt>& indices) {
//    Vector3d min{Math::Constants<Double>::inf()}, max{-Math::Constants<Double>::inf()};
//    for(auto const& p : vertices) {
//        for(int i = 0; i < 3; ++i) {
//            min[i] = Math::min(min[i], p[i]);
//            max[i] = Math::max(max[i], p[i]);
//        }
//    }
//    auto scale = (max - min).lengthInverted();
//    for(auto& v : vertices) {
//        v = (v - min)*scale;
//    }
//
//    ArrayView<Vector3d> view{vertices, vertices.size()};
//    auto afterRemoving = MeshTools::removeDuplicatesFuzzyIndexedInPlace(indices,
//                                                                        arrayCast<2, double>(view));
//    arrayResize(vertices, afterRemoving);
//}

//void normalizeMesh(Trade::MeshData& meshData) {
//    Array<char> is(NoInit, sizeof(UnsignedInt)*meshData.indexCount()),
//            vs(NoInit, sizeof(Vector3)*meshData.vertexCount());
//    auto indices = arrayCast<UnsignedInt>(is);
//    Utility::copy(meshData.indices<UnsignedInt>(), indices);
//
//    auto vertices = arrayCast<Vector3>(vs);
//    auto points = meshData.attribute<Vector3>(Trade::MeshAttribute::Position);
//
//    Vector3 min{Math::Constants<Float>::inf()}, max{-Math::Constants<Float>::inf()};
//    for(auto const& p : points) {
//        for(int i = 0; i < 3; ++i) {
//            min[i] = Math::min(min[i], p[i]);
//            max[i] = Math::max(max[i], p[i]);
//        }
//    }
//    auto scale = (max - min).lengthInverted();
//    auto mid = (max - min)*0.5f;
//    for(int i = 0; i < vertices.size(); ++i) {
//        vertices[i] = (points[i])*scale;
//    }
//
//    ArrayView<Vector3> view{vertices, vertices.size()};
//    auto afterRemoving = MeshTools::removeDuplicatesFuzzyIndexedInPlace(indices,
//                                                                        arrayCast<2, double>(view));
//
//    arrayResize(vs, afterRemoving*sizeof(Vector3));
//    vertices = arrayCast<Vector3>(vs);
//
//    Trade::MeshIndexData indexData{indices};
//    Trade::MeshAttributeData vertexData{Trade::MeshAttribute::Position, vertices};
//
//    meshData = Trade::MeshData(MeshPrimitive::Triangles, std::move(is), indexData, std::move(vs), {vertexData});
//}

void Viewer::drawBrushOptions() {
    if(ImGui::TreeNode("Brush")) {
        constexpr int step = 1;
        constexpr double stepDist = 0.01;
        constexpr double min = 0.f, max = 1.f;
        ImGui::SliderScalar("Recursive Phase Filter Factor", ImGuiDataType_Double, &recursiveFilterFactor, &min, &max,
                            "%.3f", 2.0f);
        ImGui::SliderScalar("Distance Step", ImGuiDataType_Double, &distStep, &min, &max, "%.5f", 2.0f);
        ImGui::InputScalar("Maximal Distance", ImGuiDataType_Double, &maxDist, &stepDist, nullptr, "%.3f");
        constexpr double lower = -1.f, upper = 1.f;
        ImGui::SliderScalar("Phase", ImGuiDataType_Double, &phase, &lower, &upper, "%.2f", 1.0f);

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

bool Viewer::drawFunctionals(Array<Functional>& functionals, size_t& id) {
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

        if(ImGui::Checkbox("Disable", &f.disable)) {
            proxy.redraw();
        }

        f.drawImGuiOptions(proxy);

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
            arrayAppend(problem.objectives, makeFunctional(currentType));

        ImGui::SameLine();

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
                proxy.setCallbacks([this](Viewer*) {
                    auto parameters = currentNode.phasefield();
                    Array<double> gradP(parameters.size());
                    double cost = 0;
                    problem.nodeToOptimize = currentNode;
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

        if(dirichletScaling.refCount() > 1 || doubleWellScaling.refCount() > 1 || connectednessScaling.refCount() > 1) {
            constexpr Double minEps = 0.f, maxEps = 0.1;
            if(ImGui::DragScalar("epsilon", ImGuiDataType_Double, &epsilon, .0001f, &minEps, &maxEps, "%f", 1))
                setScalingFactors();
        }

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

        ImGui::InputScalar("Max Depth", ImGuiDataType_U64, &maximumDepth, &step, nullptr, "%u");
        ImGui::Checkbox("Hierarchical Optimization", &hierarchicalOptimization);

        if(ImGui::Button("Optimize") && !problem.objectives.empty() && !isOptimizing) {
            isOptimizing = true;
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

void Viewer::runOptimization(UniqueFunction<bool()>&& cb){
    auto abortCb = [this, cb = std::move(cb)] (Solver::IterationSummary const&) mutable -> Solver::Status::Value {
        proxy.redraw();
        return cb() && isOptimizing ? Solver::Status::CONTINUE : Solver::Status::USER_ABORTED;
    };
    auto& callbacks = options.callbacks;
    arrayAppend(callbacks, InPlaceInit, std::move(abortCb));

    if(hierarchicalOptimization) {
        tree.root().initializePhasefieldFromParent();

        double largeScaleEpsilon = 0.1;
        for(size_t d = 0; d <= maximumDepth; ++d) {
            tree.computeWeightsOfAncestorsOfLevel(d);

            for(Node node : tree.nodesOnLevel(d)) {
                double area = 0.5*node.integrateWeight(mesh);
                for(Functional& f : problem.objectives) {
                    if(f.functionalType == FunctionalType::AreaRegularizer)
                        reinterpret_cast<AreaRegularizer*>(f.erased)->targetArea = area;
                }

                problem.nodeToOptimize = node;
                //double oldEpsilon = epsilon;
                //epsilon = largeScaleEpsilon;
                //setScalingFactors();
                //Solver::solve(options, problem, node.phasefield(), nullptr);
                //epsilon = oldEpsilon;
                //setScalingFactors();
                Solver::solve(options, problem, node.phasefield(), nullptr);
            }

            largeScaleEpsilon *= 0.5;
            if(d != maximumDepth) {
                for(Node node : tree.nodesOnLevel(d))
                    node.splitAndInitialize(&currentNode);
            }
        }
    } else {
        tree.computeWeightsOfAncestorsOfLevel(currentNode.depth());
        double area = 0.5*currentNode.integrateWeight(mesh);
        Debug{} << "Target area" << area;
        for(Functional& f : problem.objectives) {
            if(f.functionalType == FunctionalType::AreaRegularizer) {
                reinterpret_cast<AreaRegularizer*>(f.erased)->targetArea = area;
            }
        }
        problem.nodeToOptimize = currentNode;
        Solver::solve(options, problem, currentNode.phasefield(), nullptr);
    }

    arrayResize(callbacks, callbacks.size() - 1); /* pop the last callback off */
}

Functional Viewer::makeFunctional(FunctionalType::Value type) {
    switch(type) {
        case FunctionalType::AreaRegularizer: {
            Functional f = AreaRegularizer{mesh};
            f.loss = QuadraticLoss{};
            return f;
        }
        case FunctionalType::DirichletEnergy: {
            Functional f = DirichletEnergy{mesh};
            f.scaling = dirichletScaling;
            return f;
        }
        case FunctionalType::DoubleWellPotential : {
            Functional f = DoubleWellPotential{mesh};
            f.scaling = doubleWellScaling;
            return f;
        }
        case FunctionalType::ConnectednessConstraint : {
            Functional f = ConnectednessConstraint{mesh};
            f.scaling = connectednessScaling;
            return f;
        }
        case FunctionalType::WeakYamabe : return WeakYamabe{mesh};

        default: return Functional{};
    }
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
    *dirichletScaling = epsilon/2.;
    *doubleWellScaling = 1./epsilon;
    *connectednessScaling = 1./(epsilon*epsilon*epsilon);
}

Vector3 Viewer::unproject(Vector2i const& windowPosition) {
    auto fbSize = framebufferSize();
    auto wSize = windowSize();

    const Vector2i position = windowPosition*Vector2{fbSize}/Vector2{wSize};
    const Vector2i fbPosition{position.x(), Mg::GL::defaultFramebuffer.viewport().sizeY() - position.y() - 1};

    Float depth;
    Mg::Image2D data = GL::defaultFramebuffer.read(
            Range2Di::fromSize(fbPosition, Vector2i{1}).padded(Vector2i{2}),
            {GL::PixelFormat::DepthComponent, Mg::GL::PixelType::Float});
    depth = Math::min<Float>(arrayCast<const Float>(data.data()));

    const Vector2i viewSize = wSize;
    const Vector2i viewPosition{windowPosition.x(), viewSize.y() - windowPosition.y() - 1};
    const Vector3 in{2*Vector2{viewPosition}/Vector2{viewSize} - Vector2{1.0f}, depth*2.0f - 1.0f};

    //get global coordinates
    return (arcBall->transformationMatrix()*projection.inverted()).transformPoint(in);
}

void Viewer::drawVisualizationOptions() {
    if(ImGui::TreeNode("Visualization Options")) {

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

        static bool drawLogScaling = false;
        if(ImGui::Checkbox("Compute Optimal Conformal Scaling", &drawLogScaling)) {
            if(drawLogScaling) {
                proxy.setCallbacks(
                        [this](Viewer*){
                            VertexData<double> solution;
                            solveDiffuseYamabeEquation(mesh, currentNode.phasefield(), currentNode.temporary(), solution);
                            double min, max, w;
                            std::tie(min, max) = Math::minmax(ArrayView<double>(solution));
                            if(max - min < 1e-10)
                                w = 1;
                            else
                                w = max - min;
                            proxy.drawValues(solution, [=](double value){ return (value - min)/w; });
                            },
                        [this]{ drawLogScaling = false; });
            } else proxy.setDefaultCallback();
            proxy.redraw();
        }

        char current[100];
        sprintf(current, "%zu", currentNode.idx);
        if(ImGui::BeginCombo("##currentnode", current)) {
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


        if(ImGui::Button("Split Leaf nodes")) {
            for(Node node : tree.leafs()) {
                node.addRightChild(&currentNode);
                node.addLeftChild(&currentNode);
            }
            for(Node node : tree.leafs())
                node.initializePhasefieldFromParent();
            proxy.redraw();
        }

        if(ImGui::Button("Initialize Current Node")) {
            currentNode.initializePhasefieldFromParent();
            proxy.redraw();
        }

        ImGui::SameLine();

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

        if(ImGui::Button("Colorize With Selected Phase")) {
            for(double& x : currentNode.phasefield()) x = 0;
            proxy.redraw();
        }

        ImGui::TreePop();
    }
}

void Viewer::drawIO() {
    if(ImGui::TreeNode("IO")) {

        const char* primitives[] = {
            "capsule3DSolid",
            "circle3DSolid",
            "coneSolid",
            "cylinderSolid",
            "grid3DSolid",
            "icosphereSolid",
            "planeSolid",
            "squareSolid",
            "uvSphereSolid"
        };

        static const char* currentPrimitive = primitives[5];
        if(ImGui::BeginCombo("##primitive", currentPrimitive)) {
            for(const char* prim : primitives) {
                bool isSelected = prim == currentPrimitive;
                if(ImGui::Selectable(prim, isSelected)) {
                    currentPrimitive = prim;
                }
                if(isSelected)
                    ImGui::SetItemDefaultFocus();
            }
            ImGui::EndCombo();
        }

        static char path[50] = "/home/janos/meshes/experiments";
        static char postfix[25] = "spot";

        ImGui::PushItemWidth(150);

        ImGui::BeginGroup();
        ImGui::Text("Import/Export to");
        ImGui::InputText("##Export to", path, sizeof(path));
        ImGui::EndGroup();

        ImGui::SameLine();

        ImGui::BeginGroup();
        ImGui::Text("Mesh Name");
        ImGui::InputText("##Mesh Name", postfix, sizeof(postfix));
        ImGui::EndGroup();

        ImGui::PopItemWidth();
        ImGui::Dummy(ImVec2(0.0f, 20.0f));

        if(ImGui::Button("Export Scene")) {
            if(!Cr::Utility::Directory::exists(path))
                Cr::Utility::Directory::mkpath(path);

            Array<char> data;
            tree.serialize(data);

            if(strlen(postfix)) {
                std::string treePath = Cr::Utility::Directory::join(path, std::string{postfix} + ".bin");
                FILE* fp = fopen(treePath.c_str(), "w");
                fwrite(data.data(), sizeof(char), data.size(), fp);
                fclose(fp);

                std::string meshPath = Cr::Utility::Directory::join(path, std::string{postfix} + ".ply");
                saveMesh(meshPath.c_str(), mesh);
            }
        }

        ImGui::SameLine();
        if(ImGui::Button("Load Scene")) {
            loadScene(path, postfix);
        }

        ImGui::SameLine();
        if(ImGui::Button("Load Primitive")) {
            loadScene("deadbeef", currentPrimitive);
        }

        ImGui::Checkbox("Animate", &animate);

        static char recordingPath[100] = "/home/janos/test.mp4";
        ImGui::InputText("Video Path", recordingPath, 100);

        const auto green = ImVec4(0.56f, 0.83f, 0.26f, 1.0f);
        const auto red = ImVec4(0.83f, 0.56f, 0.26f, 1.0f);

        ImGui::PushStyleColor(ImGuiCol{}, recording ? red : green);
        if(ImGui::Button("Start Recording")) {
            if(!recording) {
                videoSaver.startRecording(recordingPath, framebufferSize());
                recording = true;
            }
        }
        ImGui::PopStyleColor();
        ImGui::SameLine();
        ImGui::PushStyleColor(ImGuiCol{}, recording ? green : red);
        if(ImGui::Button("Stop Recording")) {
            if(recording) {
                videoSaver.endRecording();
                recording = false;
            }
        }
        ImGui::PopStyleColor();

        ImGui::TreePop();
    }
}

void Viewer::viewportEvent(ViewportEvent& event) {
    Mg::GL::defaultFramebuffer.setViewport({{}, event.framebufferSize()});
    arcBall->reshape(event.windowSize());

    imgui.relayout(Vector2{event.windowSize()}/event.dpiScaling(),
                   event.windowSize(), event.framebufferSize());
}

void Viewer::keyPressEvent(KeyEvent& event) {
    if(imgui.handleKeyPressEvent(event)) {
        event.setAccepted();
        return;
    }

    switch(event.key()) {
        case KeyEvent::Key::L:
            if(arcBall->lagging() > 0.0f) {
                Debug{} << "Lagging disabled";
                arcBall->setLagging(0.0f);
            } else {
                Debug{} << "Lagging enabled";
                arcBall->setLagging(0.85f);
            }
            break;
        case KeyEvent::Key::R:
            arcBall->reset();
            break;
        case KeyEvent::Key::LeftCtrl :
            brushingModeEnabled = true;
            Mg::GL::defaultFramebuffer.mapForRead(Mg::GL::DefaultFramebuffer::ReadAttachment::Front);
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
        Mg::GL::defaultFramebuffer.mapForRead(Mg::GL::DefaultFramebuffer::ReadAttachment::None);
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

void Viewer::startBrushing(Vector3 const& point) {
    fastMarchingMethod.update();

    arrayResize(distances, 0);
    auto [idx, _] = kdtree.nearestNeighbor(point);
    fastMarchingMethod.setSource(Vertex{size_t(idx), &mesh});
    targetDist = 0.;
    brushing = true;
}

void Viewer::mousePressEvent(MouseEvent& event) {
    if(imgui.handleMousePressEvent(event)) {
        event.setAccepted();
        return;
    }

    if(brushingModeEnabled) {
        startBrushing(unproject(event.position()));
        event.setAccepted();
        return;
    }

    if(event.button() == MouseEvent::Button::Middle) {
        trackingMouse = true;
        ///* Enable mouse capture so the mouse can drag outside of the window */
        ///** @todo replace once https://github.com/mosra/magnum/pull/419 is in */
        SDL_CaptureMouse(SDL_TRUE);

        arcBall->initTransformation(event.position());

        event.setAccepted();
        redraw(); /* camera has changed, redraw! */
    }
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

    if(event.button() == MouseEvent::Button::Middle) {
        /* Disable mouse capture again */
        /** @todo replace once https://github.com/mosra/magnum/pull/419 is in */
        if(trackingMouse) {
            SDL_CaptureMouse(SDL_FALSE);
            trackingMouse = false;
            event.setAccepted();
        }
    }
}

void Viewer::mouseMoveEvent(MouseMoveEvent& event) {
    if(imgui.handleMouseMoveEvent(event)) {
        event.setAccepted();
        return;
    }

    if(brushingModeEnabled && brushing){
        startBrushing(unproject(event.position()));
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

    proxy.upload(); /* synchronize with gpu */

    if(animate) {
        Vector2i m = windowSize()/2;
        arcBall->initTransformation(m);
        arcBall->rotate(m + Vector2i{5, 0});
    }

    /* draw scene */
    bool camChanged = arcBall->updateTransformation();
    Matrix4 viewTf = arcBall->viewMatrix();


    phong.setProjectionMatrix(projection)
         .setTransformationMatrix(viewTf)
         .setNormalMatrix(viewTf.normalMatrix())
         .setLightPositions({{10,10,10}, {-10, -10, 10}})
         .setLightColors({Color4{0.5}, Color4{0.5}})
         .setSpecularColor(Color4{0.1})
         .draw(glMesh);

    if(recording) {
        Mg::Image2D image = GL::defaultFramebuffer.read({{},framebufferSize()}, {GL::PixelFormat::RGBA, Mg::GL::PixelType::UnsignedByte});
        videoSaver.appendFrame(std::move(image));
    }
    //if(recording) {
    //    Mg::Image2D image = GL::defaultFramebuffer.read({{},framebufferSize()}, {GL::PixelFormat::RGBA, Mg::GL::PixelType::UnsignedByte});
    //    videoSaver.appendFrame(std::move(image));
    //}

    /* draw ImGui stuff */
    //drawSubdivisionOptions();
    drawBrushOptions();
    drawOptimizationOptions();
    drawVisualizationOptions();
    drawIO();

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


}
