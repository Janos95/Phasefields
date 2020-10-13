//
// Created by janos on 08.11.19.
//

#include "Viewer.h"
#include "ImGuiWidgets.h"
#include "Enums.h"
#include "Tag.h"
#include "ModicaMortola.h"
#include "ConnectednessConstraint.h"
#include "DiffuseYamabe.h"
#include "C1Functions.h"

#include <ScopedTimer/ScopedTimer.h>

#include <Corrade/Utility/Algorithms.h>
#include <Corrade/Utility/Configuration.h>
#include <Corrade/Utility/Resource.h>
#include <Corrade/Containers/StaticArray.h>
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

#include <Magnum/Math/Matrix4.h>
#include <Magnum/Math/FunctionsBatch.h>

#include <Magnum/Shaders/Phong.h>

#include <Magnum/Trade/AbstractSceneConverter.h>
#include <Magnum/Trade/MeshData.h>
#include <Magnum/Trade/AbstractImporter.h>
#include <Magnum/Image.h>
#include <Magnum/ImageView.h>
#include <Magnum/PixelFormat.h>
#include <Magnum/DebugTools/ColorMap.h>

#include <Magnum/MeshTools/Compile.h>
#include <Magnum/MeshTools/RemoveDuplicates.h>

#include <Magnum/ImGuiIntegration/Context.hpp>
#include <MagnumPlugins/PrimitiveImporter/PrimitiveImporter.h>

#include <implot.h>

namespace Phasefield {

using namespace Mg::Math::Literals;

namespace {

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
               .setStorage(1, Mg::GL::TextureFormat::RGBA, size) // or SRGB8
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


Viewer::Viewer(Arguments const& arguments) :
        Mg::Platform::Application{arguments, Mg::NoCreate},
        fastMarchingMethod(mesh),
        tree{mesh},
        bvh{mesh},
        proxy(*this),
        problem(tree)
        //problem(tree),
        //proxy(*this),
        //segmentationTag(getTag()),
        //phasefieldTag(getTag())
{
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
            .setWindowFlags(Configuration::WindowFlag::Resizable)
            .setSize(conf.size(), dpiScaling);
        GLConfiguration glConf;
        glConf.setSampleCount(dpiScaling.max() < 2.0f ? 8 : 2);
#ifdef MAGNUM_TARGET_WEBGL
        /* Needed to ensure the canvas depth buffer is always Depth24Stencil8,
           stencil size is 0 by default, some browser enable stencil for that
           (Chrome) and some don't (Firefox) and thus our texture format for
           blitting might not always match. */
        glConf.setDepthBufferSize(24)
            .setStencilBufferSize(8);
#endif
        if(!tryCreate(conf, glConf))
            create(conf, glConf.setSampleCount(0));
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

        phongColorMap = Phong{Phong::Flag::DiffuseTexture, 2};
        phongVertexColors = Phong{Phong::Flag::VertexColor, 2};

        /* On WebGL we can't just read depth. The only possibility to read depth is
        to use a depth texture and read it from a shader, then reinterpret as
        color and write to a RGBA texture which can finally be read back using
        glReadPixels(). However, with a depth texture we can't use multisampling
        so I'm instead blitting the depth from the default framebuffer to
        another framebuffer with an attached depth texture and then processing
        that texture with a custom shader to reinterpret the depth as RGBA
        values, packing 8 bit of the depth into each channel. That's finally
        read back to the client. */
#ifdef MAGNUM_TARGET_WEBGL
        depth = GL::Texture2D{};
        depth.setMinificationFilter(GL::SamplerFilter::Nearest)
            .setMagnificationFilter(GL::SamplerFilter::Nearest)
            .setWrapping(GL::SamplerWrapping::ClampToEdge)
            /* The format is set to combined depth/stencil in hope it will match
               the browser depth/stencil format, requested in the GLConfiguration
               above. If it won't, the blit() won't work properly. */
            .setStorage(1, GL::TextureFormat::Depth24Stencil8, framebufferSize());
        depthFramebuffer = GL::Framebuffer{{{}, framebufferSize()}};
        depthFramebuffer.attachTexture(GL::Framebuffer::BufferAttachment::Depth, depth, 0);

        reinterpretDepth = GL::Renderbuffer{};
        reinterpretDepth.setStorage(GL::RenderbufferFormat::RGBA8, framebufferSize());
        reinterpretFramebuffer = GL::Framebuffer{{{}, framebufferSize()}};
        reinterpretFramebuffer.attachRenderbuffer(GL::Framebuffer::ColorAttachment{0}, reinterpretDepth);
        reinterpretShader = DepthReinterpretShader{};
        fullscreenTriangle = GL::Mesh{};
        fullscreenTriangle.setCount(3);
#endif

        colorMapTextures = makeColorMapTextures();

    }

    /* setup the problem */
    {
        arrayAppend(problem.objectives, makeFunctional(FunctionalType::DoubleWellPotential));
        arrayAppend(problem.objectives, makeFunctional(FunctionalType::DirichletEnergy));
        //arrayAppend(problem.objectives, makeFunctional(FunctionalType::DiffuseYamabe));
        //arrayAppend(problem.objectives, makeFunctional(FunctionalType::AreaRegularizer));
        //arrayAppend(problem.objectives, makeFunctional(FunctionalType::ConnectednessConstraint));

        size_t objectiveCount = problem.objectives.size();
        arrayResize(show, DirectInit, objectiveCount, true);
        arrayResize(data, objectiveCount);
    }

    /* try to load tree from disk, otherwise fallback to empty tree with three nodes */
    {

        Mg::PluginManager::Manager<Mg::Trade::AbstractImporter> manager;
        primitiveImporter = new Mg::Trade::PrimitiveImporter{manager, "PrimitiveImporter"};
        primitiveImporter->openData("deadbeef");

        //primitiveImporter = manager.instantiate("PrimitiveImporter");

        ScopedTimer t{"Loading scene from disk", true};
        loadScene("/home/janos/meshes/testing", "hill");

    }

    /* Setup ImGui, ImPlot, load a better font */
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

        ImPlot::CreateContext(); /* init implot context */

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
#ifndef CORRADE_TARGET_EMSCRIPTEN
    setSwapInterval(1);
    setMinimalLoopPeriod(16);
#endif

}

void Viewer::loadScene(const char* path, const char* postfix) {
#ifdef PHASEFIELD_WITH_ASSIMP
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

    } else { /* try a primitive */
        Cr::Utility::Configuration conf{"/home/janos/Phasefield/primitives.conf"};
        primitiveImporter->configuration() = *conf.group("configuration");

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
#endif
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

        ImGui::PushItemWidth(150);
        f.drawImGuiOptions(proxy);
        ImGui::PopItemWidth();

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
                proxy.setCallbacks([this](Node) {
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

        constexpr Double minEps = 0.f, maxEps = 0.1;
        if(ImGui::DragScalar("epsilon", ImGuiDataType_Double, &epsilon, .0001f, &minEps, &maxEps, "%f", 1))
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

    /* callbacks */

    auto abortCb = [this, cb = std::move(cb)] (Solver::IterationSummary const&) mutable -> Solver::Status::Value {
        proxy.redraw();
        return cb() && isOptimizing ? Solver::Status::CONTINUE : Solver::Status::USER_ABORTED;
    };
    auto& callbacks = options.callbacks;
    arrayAppend(callbacks, InPlaceInit, std::move(abortCb));

    for(auto& d : data) d.clear();
    arrayAppend(callbacks, InPlaceInit, [this](Solver::IterationSummary const& summary){
        auto& objs = problem.objectives;
        Node node = problem.nodeToOptimize;
        for(size_t i = 0; i < objs.size(); ++i) {
            double cost = 0;
            objs[i](node.phasefield(), node.temporary(), cost, nullptr, nullptr);
            data[i].add(t, cost);
        }
        ++t;
        return Solver::Status::CONTINUE;
    });

    if(hierarchicalOptimization) {
        tree.root().initializePhasefieldFromParent();

        for(size_t d = 0; d <= maximumDepth; ++d) {
            tree.computeWeightsOfAncestorsOfLevel(d);

            for(Node node : tree.nodesOnLevel(d)) {
                setAreaConstraint(node);
                problem.nodeToOptimize = node;
                Solver::solve(options, problem, node.phasefield(), nullptr);
            }

            if(d != maximumDepth) {
                for(Node node : tree.nodesOnLevel(d))
                    node.splitAndInitialize(&currentNode);
            }
        }
    } else {
        tree.computeWeightsOfAncestorsOfLevel(currentNode.depth());
        setAreaConstraint(currentNode);
        problem.nodeToOptimize = currentNode;
        Solver::solve(options, problem, currentNode.phasefield(), nullptr);
    }

    arrayResize(callbacks, callbacks.size() - 2); /* pop the last two callbacks off */

}

Functional Viewer::makeFunctional(FunctionalType::Value type) {
    switch(type) {
        case FunctionalType::AreaRegularizer: {
            Functional f = AreaRegularizer{mesh};
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
            Functional f = ConnectednessConstraint{mesh};
            f.scaling = &connectednessScaling;
            return f;
        }
        case FunctionalType::DiffuseYamabe : {
            DiffuseYamabe yamabe{mesh};
            yamabe.lambdaScaling = &yamabeLambdaScaling;
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
    dirichletScaling = epsilon*epsilon/2.;
    connectednessScaling = 1./(epsilon*epsilon);
    areaPenaltyScaling = 1.;
    doubleWellScaling = 1.;
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

        if(ImGui::Button("Print Information")) {
            SmootherStep chi;
            Array<double> patchAreas;
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
                arrayAppend(patchAreas, {posArea, negArea});
            }

            CORRADE_ASSERT(patchAreas.size() == tree.numLeafs*2, "weird number of patches",);
            double targetArea = 0;
            for(double x : patchAreas) targetArea += x;
            targetArea /= double(patchAreas.size());
            printf("Target Area %f\n", targetArea);
            double totalError = 0;
            for(double patchArea : patchAreas) {
                double error = Math::abs(targetArea - patchArea);
                totalError += error;
                printf("Patch Error %f, (area = %f)\n", error, patchArea);
            }
            printf("Total error %f\n", totalError);
        }

        static bool drawCurvature = false;
        if(ImGui::Checkbox("Gaussian Curvature", &drawCurvature)) {
            if(drawCurvature) {
                proxy.setCallbacks(
                        [this](Node node) {
                            mesh.requireGaussianCurvature();
                            double min = std::numeric_limits<double>::max();
                            double max = std::numeric_limits<double>::min();
                            for(Vertex v : mesh.vertices()) {
                                if(!v.onBoundary()) {
                                    min = Math::min(min, mesh.gaussianCurvature[v]);
                                    max = Math::max(max, mesh.gaussianCurvature[v]);
                                }
                            }
                            Debug{} << min;
                            Debug{} << max;
                            proxy.drawValues(mesh.gaussianCurvature, [=](double x) { return (x-min)/(max-min); });
                        },
                        [this]{ drawCurvature = false; });
            } else proxy.setDefaultCallback();
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
#ifdef PHASEFIELD_WITH_FFMPEG
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
#ifdef PHASEFIELD_WITH_FFMPEG
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

    imgui.relayout(Vector2{event.windowSize()}/event.dpiScaling(),
                   event.windowSize(), event.framebufferSize());

    /* Recreate textures and renderbuffers that depend on viewport size */
#ifdef MAGNUM_TARGET_WEBGL
    depth = GL::Texture2D{};
    depth.setMinificationFilter(GL::SamplerFilter::Nearest)
        .setMagnificationFilter(GL::SamplerFilter::Nearest)
        .setWrapping(GL::SamplerWrapping::ClampToEdge)
        .setStorage(1, GL::TextureFormat::Depth24Stencil8, event.framebufferSize());
    depthFramebuffer.attachTexture(GL::Framebuffer::BufferAttachment::Depth, depth, 0);

    reinterpretDepth = GL::Renderbuffer{};
    reinterpretDepth.setStorage(GL::RenderbufferFormat::RGBA8, event.framebufferSize());
    reinterpretFramebuffer.attachRenderbuffer(GL::Framebuffer::ColorAttachment{0}, reinterpretDepth);

    reinterpretFramebuffer.setViewport({{}, event.framebufferSize()});
#endif
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
    size_t idx = bvh.computeIntersection(origin, dir);
    if(idx == Invalid) idx = lastIntersectionIdx;
    else lastIntersectionIdx = idx;
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
        Vector3 o = unproject(event.position(), 0);
        Vector3 d = unproject(event.position(), 1) - o;
        startBrushing(o, d);
        event.setAccepted();
        return;
    }

    if(event.button() == MouseEvent::Button::Middle) {
        trackingMouse = true;
        ///* Enable mouse capture so the mouse can drag outside of the window */
        ///** @todo replace once https://github.com/mosra/magnum/pull/419 is in */

#ifndef CORRADE_TARGET_EMSCRIPTEN
        SDL_CaptureMouse(SDL_TRUE);
#endif

        arcBall->initTransformation(event.position());

        event.setAccepted();
        redraw(); /* camera has changed, redraw! */
    }

    if(event.button() == MouseEvent::Button::Right) {
        Vector3 o = unproject(event.position(), 0);
        Vector3 d = unproject(event.position(), 1) - o;

        size_t idx = bvh.computeIntersection(o, d);
        if(idx == Invalid) idx = lastIntersectionIdx;
        else lastIntersectionIdx = idx;
        Debug{} << "Phasefield value at mouse location" << currentNode.phasefield()[idx];
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
#ifndef CORRADE_TARGET_EMSCRIPTEN
            SDL_CaptureMouse(SDL_FALSE);
#endif
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

    proxy.upload(); /* synchronize with gpu */

    if(animate) {
        Vector2i m = windowSize()/2;
        arcBall->initTransformation(m);
        arcBall->rotate(m + Vector2i{5, 0});
    }

    /* draw scene */
    bool camChanged = arcBall->updateTransformation();
    Matrix4 viewTf = arcBall->viewMatrix();

    if(proxy.shaderConfig == VisualizationProxy::ShaderConfig::ColorMaps) {
        phongColorMap.setProjectionMatrix(projection)
                     .setTransformationMatrix(viewTf)
                     .setNormalMatrix(viewTf.normalMatrix())
                     .bindDiffuseTexture(colorMapTextures[colorMapIndex].second)

                     .draw(glMesh);
    } else if (proxy.shaderConfig == VisualizationProxy::ShaderConfig::VertexColors) {
        phongVertexColors.setProjectionMatrix(projection)
                         .setTransformationMatrix(viewTf)
                         .setNormalMatrix(viewTf.normalMatrix())
                         .setLightPositions({{10,10,10}, {-10, -10, 10}})
                         .setLightColors({Color4{0.5}, Color4{0.5}})
                         .setSpecularColor(Color4{0.1})
                         .draw(glMesh);
    }


    if(recording) {
#ifdef PHASEFIELD_WITH_FFMPEG
        Mg::Image2D image = GL::defaultFramebuffer.read({{},framebufferSize()}, {GL::PixelFormat::RGBA, Mg::GL::PixelType::UnsignedByte});
        videoSaver.appendFrame(std::move(image));
#endif
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

void Viewer::setAreaConstraint(Node node) {
    double area = node.integrateWeight(mesh);
    for(Functional& f : problem.objectives) {
        if(f.functionalType == FunctionalType::AreaRegularizer)
            reinterpret_cast<AreaRegularizer*>(f.erased)->totalArea = area;
    }
    for(Functional& f : problem.constraints) {
        if(f.functionalType == FunctionalType::AreaRegularizer)
            reinterpret_cast<AreaRegularizer*>(f.erased)->totalArea = area;
    }
}

Viewer::~Viewer() {
    ImPlot::DestroyContext();
}

void Viewer::drawErrorPlot() {
    bool opened;
    ImGui::Begin("ImPlot Demo", &opened, ImGuiWindowFlags_MenuBar);

    size_t objectiveCount = problem.objectives.size();

    ImGui::SameLine();

    if (ImGui::Button(paused ? "Resume" : "Pause", ImVec2(100,0)))
        paused = !paused;


    ImPlot::SetNextPlotLimitsX((double)t - 10, t, paused ? ImGuiCond_Once : ImGuiCond_Always);
    if (ImPlot::BeginPlot("##DND", nullptr, nullptr, ImVec2(-1,0), ImPlotFlags_YAxis2 | ImPlotFlags_YAxis3, ImPlotAxisFlags_NoTickLabels)) {
        for (int i = 0; i < objectiveCount; ++i) {
            const char* label = FunctionalType::to_string(problem.objectives[i].functionalType);
            if (show[i] && data[i].size() > 0) {
                ImPlot::SetPlotYAxis(0);
                ImPlot::PlotLine(label, &data[i].data[0].x(), &data[i].data[0].y(), data[i].data.size(), data[i].offset, 2 * sizeof(float));
                // allow legend labels to be dragged and dropped
                if (ImPlot::BeginLegendDragDropSource(label)) {
                    ImGui::SetDragDropPayload("DND_PLOT", &i, sizeof(int));
                    ImGui::TextUnformatted(label);
                    ImPlot::EndLegendDragDropSource();
                }
            }
        }

        ImPlot::EndPlot();
    }
    ImGui::End();
}

}
