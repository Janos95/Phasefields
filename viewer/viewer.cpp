//
// Created by janos on 08.11.19.
//

#include "viewer.hpp"
#include "subdivision.hpp"
#include "primitives.hpp"
#include "upload.hpp"
#include "optimization_context.hpp"
#include "modica_mortola.hpp"
#include "geodesic_algorithm_exact.h"
#include "connectedness_constraint.hpp"
#include "custom_widgets.hpp"

#include <scoped_timer/scoped_timer.hpp>

#include <Corrade/Utility/Algorithms.h>

#include <Magnum/GL/DefaultFramebuffer.h>
#include <Magnum/GL/Renderer.h>
#include <Magnum/GL/PixelFormat.h>
#include <Magnum/GL/TextureFormat.h>
#include <Magnum/GL/Buffer.h>
#include <Magnum/GL/DefaultFramebuffer.h>
#include <Magnum/GL/Mesh.h>

#include <Magnum/Math/Matrix4.h>
#include <Magnum/Math/FunctionsBatch.h>

#include <Magnum/Shaders/Flat.h>
#include <Magnum/Shaders/MeshVisualizer.h>
#include <Magnum/Shaders/Phong.h>
#include <Magnum/Shaders/VertexColor.h>
#include <Magnum/SceneGraph/Object.h>
#include <Magnum/ImGuiIntegration/Context.hpp>
#include <Magnum/Trade/MeshData.h>
#include <Magnum/Image.h>
#include <Magnum/PixelFormat.h>
#include <Magnum/ImageView.h>
#include <Magnum/DebugTools/ColorMap.h>


using namespace Corrade;
using namespace Magnum;

using namespace Math::Literals;

std::unordered_map<ShaderType, Cr::Containers::Pointer<Mg::GL::AbstractShaderProgram>> makeShaders() {
    std::unordered_map<ShaderType, Cr::Containers::Pointer<Mg::GL::AbstractShaderProgram>> map;

    auto Wireframe = Shaders::MeshVisualizer3D::Flag::Wireframe;
    auto PrimitiveColored = Shaders::MeshVisualizer3D::Flag::PrimitiveId;

    map.emplace(ShaderType::PhongDiffuse, new Shaders::Phong(Shaders::Phong::Flag::DiffuseTexture));
    map.emplace(ShaderType::MeshVisualizer, new Shaders::MeshVisualizer3D(Wireframe));
    map.emplace(ShaderType::MeshVisualizerPrimitiveId, new Shaders::MeshVisualizer3D(PrimitiveColored));
    map.emplace(ShaderType::FlatTextured, new Shaders::Flat3D(Shaders::Flat3D::Flag::Textured));
    map.emplace(ShaderType::VertexColor, new Shaders::VertexColor3D());

    return map;
}


std::unordered_map<ColorMapType, Mg::GL::Texture2D> makeColorMapTextures(){
    std::unordered_map<ColorMapType, Mg::GL::Texture2D> map;
    using L = std::initializer_list<std::pair<ColorMapType, Containers::StaticArrayView<256, const Vector3ub>>>;
    for(auto&& [type, colorMap] : L{
            {ColorMapType::Turbo, Magnum::DebugTools::ColorMap::turbo()},
            {ColorMapType::Magma, Magnum::DebugTools::ColorMap::magma()},
            {ColorMapType::Plasma, Magnum::DebugTools::ColorMap::plasma()},
            {ColorMapType::Inferno, Magnum::DebugTools::ColorMap::inferno()},
            {ColorMapType::Viridis, Magnum::DebugTools::ColorMap::viridis()}
    })
    {
        const Magnum::Vector2i size{Magnum::Int(colorMap.size()), 1};
        GL::Texture2D texture;
        texture.setMinificationFilter(Magnum::SamplerFilter::Linear)
                .setMagnificationFilter(Magnum::SamplerFilter::Linear)
                .setWrapping(Magnum::SamplerWrapping::ClampToEdge) // or Repeat
                .setStorage(1, Magnum::GL::TextureFormat::RGB8, size) // or SRGB8
                .setSubImage(0, {}, ImageView2D{Magnum::PixelFormat::RGB8Srgb, size, colorMap});
        map.emplace(type, std::move(texture));
    }
    return map;
}

solver::Status OptimizationCallback::operator ()(solver::IterationSummary const&){
    if(!optimize) return solver::Status::ABORT;
    return solver::Status::CONTINUE;
}


Viewer::Viewer(int argc, char** argv):
    Platform::Application{{argc,argv},NoCreate},
    dirichletScaling(1.), doubleWellScaling(1.), connectednessScaling(1.),
    optimizationCallback{optimizing}
{
    {
        Containers::arrayAppend(options.callbacks, Containers::InPlaceInit, optimizationCallback);
        problem.meshData = &meshData;
        problem.update = &update;
        problem.mutex = &mutex;
    }
    /* Setup window */
    {
        const Vector2 dpiScaling = this->dpiScaling({});
        Configuration conf;
        conf.setTitle("Viewer")
                .setSize(conf.size(), dpiScaling)
                .setWindowFlags(Configuration::WindowFlag::Resizable);
        GLConfiguration glConf;
        glConf.setSampleCount(dpiScaling.max() < 2.0f ? 8 : 2);
        if(!tryCreate(conf, glConf)) {
            create(conf, glConf.setSampleCount(0));
        }
    }

    /* setup shaders and color map textures */
    {
        auto m = Primitives::cylinderSolid(10, 10, 10.f);
        auto ps = m.mutableAttribute<Vector3>(Trade::MeshAttribute::Position);
        mesh = GL::Mesh{};
        vertexBuffer = GL::Buffer{};
        indexBuffer = GL::Buffer{};
        colorMapTextures = makeColorMapTextures();
        shaders = makeShaders();
    }

    /* Setup ImGui, load a better font */
    {
        ImGui::CreateContext();
        ImGui::StyleColorsDark();

        ImFontConfig fontConfig;
        fontConfig.FontDataOwnedByAtlas = false;
        const Vector2 size = Vector2{windowSize()}/dpiScaling();
        Utility::Resource rs{"fonts"};
        Containers::ArrayView<const char> font = rs.getRaw("SourceSansPro-Regular.ttf");
        ImGui::GetIO().Fonts->AddFontFromMemoryTTF(
                const_cast<char*>(font.data()), Int(font.size()),
                20.0f*framebufferSize().x()/size.x(), &fontConfig);

        imgui = ImGuiIntegration::Context{*ImGui::GetCurrentContext(),
                                            Vector2{windowSize()}/dpiScaling(), windowSize(), framebufferSize()};

        /* Setup proper blending to be used by ImGui */
        GL::Renderer::setBlendEquation(GL::Renderer::BlendEquation::Add,
                                       GL::Renderer::BlendEquation::Add);
        GL::Renderer::setBlendFunction(GL::Renderer::BlendFunction::SourceAlpha,
                                       GL::Renderer::BlendFunction::OneMinusSourceAlpha);

    }

    /* Setup the arcball after the camera objects */
    {
        const Vector3 eye = Vector3::zAxis(-10.0f);
        const Vector3 center{};
        const Vector3 up = Vector3::yAxis();
        camera.emplace(scene, eye, center, up, 45.0_degf,
                       windowSize(), framebufferSize());
    }

    object = new Object3D(&scene);

    GL::Renderer::enable(GL::Renderer::Feature::DepthTest);
    GL::Renderer::enable(GL::Renderer::Feature::FaceCulling);

    /* Start the timer, loop at 60 Hz max */
    setSwapInterval(1);
    setMinimalLoopPeriod(16);
}



void Viewer::drawSubdivisionOptions() {
    if (ImGui::TreeNode("Subdivisions")) {
        ImGui::Text("Currently we have %d vertices and %d faces", (int)vertices.size(), (int)indices.size() / 3);
        constexpr int step = 1;
        ImGui::InputScalar("Number of Sibdivision (wrt. orignal mesh)", ImGuiDataType_U32, &numSubdivisions, &step, nullptr, "%d");
        if (ImGui::Button("do subdivision")) {
            subdivide(numSubdivisions, original, phasefield, meshData);
            updateInternalDataStructures();
            upload(mesh, vertexBuffer, indexBuffer, meshData);
            updateFunctionals(problem.functionals);
            updateFunctionals(problem.constraints);
            makeDrawableCurrent(drawableType);
        }
        ImGui::TreePop();
    }
}

void Viewer::makeDrawableCurrent(DrawableType type) {
    if(!object) return;
    object->features().clear();
    drawableType = type;
    switch(type){
        case DrawableType::FaceColored : {
            auto d = new FaceColorDrawable(*object, mesh, *shaders[ShaderType::MeshVisualizerPrimitiveId], &drawableGroup);
            d->texture = faceTexture.get();
            d->offset = 0; d->scale = 3.f/static_cast<Float>(indices.size());
            drawable = d;
            break;
        }
        case DrawableType::PhongDiffuse : {
            auto d = new PhongDiffuseDrawable(*object, mesh, *shaders[ShaderType::PhongDiffuse], &drawableGroup);
            d->texture = &colorMapTextures[colorMapType];
            drawable = d;
            break;
        }
        case DrawableType::FlatTextured : {
            auto d = new FlatDrawable(*object, mesh, *shaders[ShaderType::FlatTextured], &drawableGroup);
            d->texture = &colorMapTextures[colorMapType];
            drawable = d;
            break;
        }
        case DrawableType::MeshVisualizer : {
            drawable = new MeshVisualizerDrawable(*object, mesh, *shaders[ShaderType::MeshVisualizer], &drawableGroup);
            break;
        }
    }
}

void Viewer::drawPrimitiveOptions() {
    if(handlePrimitive(original, expression)){
        meshData = preprocess(original, CompileFlag::GenerateSmoothNormals|CompileFlag::AddTextureCoordinates);
        updateInternalDataStructures();
        upload(mesh, vertexBuffer, indexBuffer, meshData);
        updateFunctionals(problem.functionals);
        updateFunctionals(problem.constraints);
        makeDrawableCurrent(drawableType); //for first time
    }
}

void Viewer::drawBrushOptions() {
    if (ImGui::TreeNode("Brush"))
    {
        constexpr int step = 1;
        constexpr double stepDist = 0.01;
        constexpr double min = 0.f, max = 1.f;
        ImGui::SliderScalar("Recursive Phase Filter Factor", ImGuiDataType_Double, &recursiveFilterFactor, &min, &max, "%.3f", 2.0f);
        ImGui::SliderScalar("Distance Step", ImGuiDataType_Double, &distStep, &min, &max, "%.5f", 2.0f);
        ImGui::InputScalar("Maximal Distance", ImGuiDataType_Double, &maxDist, &stepDist, nullptr, "%.3f");
        constexpr double lower = -1.f, upper = 1.f;
        ImGui::SliderScalar("Phase", ImGuiDataType_Double, &phase, &lower, &upper, "%.2f", 1.0f);

        const auto colBrushing = ImVec4(0.56f, 0.83f, 0.26f, 1.0f);
        const auto colNotBrushing = ImVec4(0.85f, 0.85f, 0.85f, 1.0f);
        ImGui::TextColored(brushing ? colBrushing : colNotBrushing, "Press Left Control To Enable Brushing");
        ImGui::SameLine();
        ImVec2 p = ImGui::GetCursorScreenPos();
        ImDrawList* draw_list = ImGui::GetWindowDrawList();
        float height = ImGui::GetFrameHeight();
        float width = height * 1.55f;
        ImGui::InvisibleButton("brush modus", ImVec2(width, height));
        draw_list->AddRectFilled(p, ImVec2(p.x + width, p.y + height), brushing ? ImGui::GetColorU32(colBrushing) : ImGui::GetColorU32(colNotBrushing), height * 0.1f);
        ImGui::TreePop();
    }
}

void Viewer::makeExclusiveVisualizer(Functional* funcPtr) {
    if(exclusiveVisualizer && exclusiveVisualizer != funcPtr->metaData.get()) {
        std::lock_guard l(mutex);
        exclusiveVisualizer->flags = {};
    }
    exclusiveVisualizer = funcPtr->metaData.get();
    problem.flags = {};
}

bool Viewer::drawGradientMetaData(
        GradientMetaData& meta,
        bool& makeExclusive,
        bool& evaluateProblem){
    bool drawGrad = bool(meta.flags & VisualizationFlag::Gradient);
    if(ImGui::Checkbox("Gradient", &drawGrad) ){
        std::lock_guard l(mutex);
        if(drawGrad){
            makeDrawableCurrent(DrawableType::PhongDiffuse);
            meta.flags = VisualizationFlag::Gradient;
            makeExclusive = true;
            evaluateProblem = true;
        } else {
            meta.flags = {};
        };
    }
}

/**
 * returns true if an event triggered an exclusive visualizations options.
 * Also note that we can safely read from meta.flags since the optimization thread
 * only reads from those as well (even taking the lock in case we modify)
 */
void Viewer::drawConnectednessConstraintOptions(
        ConnectednessMetaData<Double>& meta,
        bool& makeExclusive,
        bool& evaluateProblem){
    double a = meta.a, b = meta.b;
    if(dragDoubleRange2("Positive Interval", &a, &b, 0.01f, -1.f, 1.f, "Min: %.2f", "Max: %.2f",1.f)){
        std::lock_guard l(mutex);
        meta.a = a; meta.b = b;
    }

    ImGui::BeginGroup();
    if(ImGui::Checkbox("Visualize Shortest Paths", &meta.paths->drawPaths) ){
        std::lock_guard l(mutex);
        if(meta.paths->drawPaths){
            evaluateProblem = true;
            meta.generateLineStrips = true;
        } else {
            meta.generateLineStrips = false;
        }
    }

    bool visGeodesicWeights = static_cast<bool>(meta.flags & VisualizationFlag::GeodesicWeights);
    if(ImGui::Checkbox("Geodesic Weights", &visGeodesicWeights)){
        std::lock_guard l(mutex);
        if(visGeodesicWeights) {
            makeDrawableCurrent(DrawableType::FaceColored);
            makeExclusive = true;
            meta.flags = VisualizationFlag::GeodesicWeights;
            evaluateProblem = true;
        }
        else meta.flags &= ~VisualizationFlag::GeodesicWeights;
    }
    ImGui::EndGroup();
    ImGui::SameLine();
    ImGui::BeginGroup();
    bool visComponents = static_cast<bool>(meta.flags & VisualizationFlag::ConnectedComponents);
    if(ImGui::Checkbox("Connected Components", &visComponents)){
        std::lock_guard l(mutex);
        if(visComponents) {
            makeDrawableCurrent(DrawableType::FaceColored);
            makeExclusive = true;
            meta.flags = VisualizationFlag::ConnectedComponents;
            evaluateProblem = true;
        }
        else meta.flags &= ~VisualizationFlag::ConnectedComponents;
    }
    bool visGradient = static_cast<bool>(meta.flags & VisualizationFlag::Gradient);
    if(ImGui::Checkbox("Gradient", &visGradient)){
        std::lock_guard l(mutex);
        if(visGradient) {
            makeDrawableCurrent(DrawableType::PhongDiffuse);
            makeExclusive = true;
            meta.flags = VisualizationFlag::Gradient;
            evaluateProblem = true;
        }
        else meta.flags &= ~VisualizationFlag::Gradient;
    }
    ImGui::EndGroup();
}

void Viewer::drawOptimizationContext() {
    if (ImGui::TreeNode("Optimization Options"))
    {
        static auto map = makeComboMapFunctionals();
        static auto cur = map.end();
        if(ImGui::BeginCombo("##combo", cur != map.end() ? cur->name.c_str() : nullptr)){
            for (auto it = map.begin(); it < map.end(); ++it){
                bool isSelected = (cur == it);
                if (ImGui::Selectable(it->name.c_str(), isSelected)){
                    cur = it;
                }
                if (isSelected)
                    ImGui::SetItemDefaultFocus();
            }
            ImGui::EndCombo();
        }

        if(ImGui::Button("Add Functional As Cost") && cur != map.end())
            Containers::arrayAppend(problem.functionals, makeFunctional(cur->type));

        ImGui::SameLine();

        if(ImGui::Button("Add Functional As Constraint") && cur != map.end())
            Containers::arrayAppend(problem.functionals, makeFunctional(cur->type));

        int nodeCount = 0;
        int toRemove = -1;
        bool evaluateProblem = false;
        auto& fs = problem.functionals;
        for(int i = 0; i < fs.size(); ++i){
            auto f = fs[i].get();
            ImGui::PushID(++nodeCount);
            ImGui::Separator();

            if(ImGui::Button("Remove")) toRemove = i;
            ImGui::SameLine();

            bool makeExclusive = false;
            constexpr double min = 0.f, max = 1.;

            switch(f->type){
                case FunctionalType::DoubleWellPotential :
                    ImGui::Text("Double Well Potential");
                    drawGradientMetaData(dynamic_cast<GradientMetaData&>(*f->metaData), makeExclusive, evaluateProblem);
                    break;
                case FunctionalType::DirichletEnergy :
                    ImGui::Text("Dirichlet Energy");
                    drawGradientMetaData(dynamic_cast<GradientMetaData&>(*f->metaData), makeExclusive, evaluateProblem);
                    break;
                case FunctionalType::Area1 : {
                    ImGui::Text("Area Regularization Quadratic");
                    drawGradientMetaData(dynamic_cast<GradientMetaData&>(*f->metaData), makeExclusive, evaluateProblem);
                    auto& farea1 = dynamic_cast<AreaRegularizer1&>(*f);
                    ImGui::SliderScalar("Area Ratio", ImGuiDataType_Double, &farea1.areaRatio, &min, &max, "%.5f", 1.0f);
                    break;
                }
                case FunctionalType::Area2 : {
                    ImGui::Text("Area Regularization Smooth Step");
                    drawGradientMetaData(dynamic_cast<GradientMetaData&>(*f->metaData), makeExclusive, evaluateProblem);
                    auto& farea2 = dynamic_cast<AreaRegularizer2&>(*f);
                    ImGui::SliderScalar("Area Ratio", ImGuiDataType_Double, &farea2.areaRatio, &min, &max, "%.5f", 1.0f);
                    break;
                }
                case FunctionalType::Connectedness :
                    ImGui::Text("Connectedness Constraint");
                    drawConnectednessConstraintOptions(
                            dynamic_cast<ConnectednessMetaData<Double>&>(*f->metaData),
                            makeExclusive,
                            evaluateProblem);
                    break;
            }

            if(makeExclusive) makeExclusiveVisualizer(f);

            drawLoss(f->metaData->loss, nodeCount);
            ImGui::PopID();
        }

        if(nodeCount) ImGui::Separator();

        if(toRemove >= 0){
            auto f = problem.functionals[toRemove].get();
            if(f->metaData.get() == exclusiveVisualizer)
                exclusiveVisualizer = nullptr;

            std::swap(fs[toRemove], fs.back());
            Containers::arrayResize(fs, fs.size() - 1);
        }

        bool visPhasefield = static_cast<bool>(problem.flags & VisualizationFlag::Phasefield);
        if(ImGui::Checkbox("Phasefield", &visPhasefield)){
            if(visPhasefield) {
                makeDrawableCurrent(DrawableType::PhongDiffuse);
                std::lock_guard l(mutex);
                problem.flags = VisualizationFlag::Phasefield;
                update = VisualizationFlag::Phasefield;
                evaluateProblem = true;
                if (exclusiveVisualizer) {
                    exclusiveVisualizer->flags = {};
                    exclusiveVisualizer = nullptr;
                }
            }
        }
        ImGui::SameLine();
        bool visGradient = static_cast<bool>(problem.flags & VisualizationFlag::Gradient);
        if(ImGui::Checkbox("Gradient", &visGradient)){
            if(visGradient) {
                makeDrawableCurrent(DrawableType::PhongDiffuse);
                std::lock_guard l(mutex);
                problem.flags = VisualizationFlag::Gradient;
                update = VisualizationFlag::Gradient;
                evaluateProblem = true;
                if (exclusiveVisualizer) {
                    exclusiveVisualizer->flags = {};
                    exclusiveVisualizer = nullptr;
                }
            }
        }

        static Double epsilon = 0.3;
        if(dirichletScaling.refCount() > 1 || doubleWellScaling.refCount() > 1 || connectednessScaling.refCount() > 1){
            constexpr Mg::Double minEps = 0.f, maxEps = 1.;
            ImGui::DragScalar("epsilon", ImGuiDataType_Double, &epsilon, .01f, &minEps, &maxEps, "%f", 2);
            *dirichletScaling = epsilon/2.;
            *doubleWellScaling = 1./epsilon;
            *connectednessScaling = 1./(epsilon * epsilon);
        }

        static std::uint32_t iterations = 100;
        constexpr static std::uint32_t step = 1;
        ImGui::InputScalar("iterations", ImGuiDataType_S32, &options.max_num_iterations, &step, nullptr, "%u");

        static std::map<solver::LineSearchDirectionType, std::string> searchDirectionMapping= {
                {solver::LineSearchDirectionType::STEEPEST_DESCENT, "Steepest Descent"},
                {solver::LineSearchDirectionType::NONLINEAR_CONJUGATE_GRADIENT, "Nonlinear Conjugate Gradient"},
                {solver::LineSearchDirectionType::LBFGS, "LBFGS"},
                {solver::LineSearchDirectionType::BFGS, "BFGS"}
        };

        if (ImGui::BeginCombo("##descent direction", searchDirectionMapping[options.line_search_direction_type].c_str())) {
            for (auto& [type, name] : searchDirectionMapping) {
                bool isSelected = (options.line_search_direction_type == type);
                if (ImGui::Selectable(name.c_str(), isSelected)){
                    options.line_search_direction_type = type;
                }
                if (isSelected)
                    ImGui::SetItemDefaultFocus();
            }
            ImGui::EndCombo();
        }
        ImGui::SameLine();
        ImGui::Text("Descent Direction");


        if(evaluateProblem && !optimizing){
            double r;
            Containers::Array<Mg::Double> grad(Containers::NoInit, phasefield.size());
            problem.evaluate(phasefield.data(), &r, grad.data());
        }

        if (ImGui::Button("Optimize") && !problem.functionals.empty())
            startOptimization();

        ImGui::SameLine();

        if (ImGui::Button("Stop"))
            stopOptimization();

        ImGui::TreePop();
    }
}

void Viewer::startOptimization() {
    stopOptimization();

    optimizing = true;
    g.run([this]{
        solver::Summary summary;
        solve(options, problem, phasefield.data(), &summary);
        optimizing = false;
        //Debug{} << summary.briefReport().c_str();
    });
}

void Viewer::stopOptimization() {
    optimizing = false;
    g.wait();
}

/* @todo would be nice to rework this at some point to handle hard deps in constructors.. */
Containers::Pointer<Functional> Viewer::makeFunctional(FunctionalType type) {
    auto ts = Containers::arrayCast<Vector3ui>(indices);
    switch (type) {
        case FunctionalType::Area1 : {
            auto p = Containers::pointer<AreaRegularizer1>(vertices, ts);
            auto meta = Containers::pointer<GradientMetaData>(meshData, update, mutex);
            meta->loss = std::move(p->metaData->loss);
            p->metaData = std::move(meta) ;
            p->type = FunctionalType::Area1; //@todo somehow this does not pick up the type from ctor
            return p;
        }
        case FunctionalType::Area2 : {
            auto p = Containers::pointer<AreaRegularizer2>(vertices, ts);
            auto meta = Containers::pointer<GradientMetaData>(meshData, update, mutex);
            meta->loss = std::move(p->metaData->loss);
            p->metaData = std::move(meta);
            p->type = FunctionalType::Area2;
            return p;
        }
        case FunctionalType::DirichletEnergy : {
            auto p = Containers::pointer<DirichletEnergy>(vertices, ts);
            auto meta = Containers::pointer<GradientMetaData>(meshData, update, mutex);
            meta->loss = std::move(p->metaData->loss);
            p->metaData = std::move(meta);
            p->metaData->scaling = dirichletScaling;
            return p;
        }
        case FunctionalType::DoubleWellPotential : {
            auto p = Containers::pointer<DoubleWellPotential>(vertices, ts);
            auto meta = Containers::pointer<GradientMetaData>(meshData, update, mutex);
            meta->loss = std::move(p->metaData->loss);
            p->metaData = std::move(meta);
            p->metaData->scaling = doubleWellScaling;
            p->type = FunctionalType::DoubleWellPotential;
            return p;
        }
        case FunctionalType::Connectedness :
            auto meta = Containers::pointer<ConnectednessMetaData<Double>>(meshData, faceColors, update, mutex);
            meta->paths = new Paths(&scene, drawableGroup);
            auto p = Containers::pointer<ConnectednessConstraint<Double>>(vertices, ts, std::move(meta));
            p->metaData->scaling = connectednessScaling;
            return p;
    }
    return nullptr;
}

void Viewer::updateFunctionals(Containers::Array<Containers::Pointer<Functional>>& functionals){
    for(auto& f : functionals){
        auto meta = std::move(f->metaData);
        f = makeFunctional(f->type);
        f->metaData = std::move(meta);
    }
}

void Viewer::paint(){
    if(targetDist < maxDist)
        targetDist += distStep;

    auto textureCoords = meshData.mutableAttribute<Vector2>(Trade::MeshAttribute::TextureCoordinates);
    for(auto [d, i] : distances) {
        if(d > targetDist) break;
        auto u = (1.f - recursiveFilterFactor) * phasefield[i] + recursiveFilterFactor * phase;
        phasefield[i] = u;
        textureCoords[i].x() = .5f * (u + 1.f);
    }
}

void Viewer::geodesicSearch() {
    geodesic::Mesh geomesh;
    //@todo this does not need to be done on every new mouse click
    //maybe check out the cgal geodesic module?
    geomesh.initialize_mesh_data(Containers::arrayCast<const Double>(vertices), indices);		//create internal mesh data structure including edges
    geodesic::GeodesicAlgorithmExact exactGeodesics(&geomesh);

    Containers::arrayResize(distances, vertices.size());
    auto it = std::min_element(vertices.begin(), vertices.end(),
                               [&](auto const &v1, auto const &v2) {
                                   return (v1 - point).dot() < (v2 - point).dot();
                               });
    auto source = it - vertices.begin();
    ScopedTimer timer("Computing geodesics", true);

    std::vector<geodesic::SurfacePoint> sources;
    sources.clear();
    sources.emplace_back(&geomesh.vertices()[source]);
    exactGeodesics.propagate(sources);
    for(unsigned i=0; i<geomesh.vertices().size(); ++i) {
        geodesic::SurfacePoint p(&geomesh.vertices()[i]);
        double distance;
        exactGeodesics.best_source(p, distance);
        distances[i].first = distance;
        distances[i].second = i;
    }
    std::sort(distances.begin(), distances.end(),[](auto& a, auto& b){ return a.first < b.first; });
}

Vector3 Viewer::unproject(Vector2i const& windowPosition) {
    auto fbSize = framebufferSize();
    auto wSize = windowSize();

    const Vector2i position = windowPosition*Vector2{fbSize}/Vector2{wSize};
    const Vector2i fbPosition{position.x(), GL::defaultFramebuffer.viewport().sizeY() - position.y() - 1};

    Float depth;
    {
        ScopedTimer timer("Reading depth from framebuffer", true);
        Image2D data = GL::defaultFramebuffer.read(
                Range2Di::fromSize(fbPosition, Vector2i{1}).padded(Vector2i{2}),
                {GL::PixelFormat::DepthComponent, GL::PixelType::Float});
        depth = Math::min<Float>(Containers::arrayCast<const Float>(data.data()));
    }

    const Vector2i viewSize = wSize;
    const Vector2i viewPosition{windowPosition.x(), viewSize.y() - windowPosition.y() - 1};
    const Vector3 in{2 * Vector2{viewPosition} / Vector2{viewSize} - Vector2{1.0f}, depth * 2.0f - 1.0f};

    //get global coordinates
    return (camera->transformationMatrix() * camera->projectionMatrix().inverted()).transformPoint(in);
}

void Viewer::updateInternalDataStructures(){
    CORRADE_ASSERT(meshData.hasAttribute(Trade::MeshAttribute::TextureCoordinates),"error", );
    auto points = meshData.attribute<Vector3>(Trade::MeshAttribute::Position);
    Containers::arrayResize(vertices, points.size());
    for (int i = 0; i < points.size(); ++i)
        vertices[i] = Vector3d(points[i]);
    indices = meshData.indicesAsArray();
    upload(mesh, vertexBuffer, indexBuffer, meshData);
    const Magnum::Vector2i size(indices.size() / 3, 1);

    faceTexture.reset(new GL::Texture2D{});
    faceTexture->setStorage(1, Magnum::GL::TextureFormat::RGB8, size)
                .setMinificationFilter(Magnum::SamplerFilter::Nearest)
                .setMagnificationFilter(Magnum::SamplerFilter::Nearest)
                .setWrapping(Magnum::SamplerWrapping::ClampToEdge);

    Containers::arrayResize(faceColors, indices.size() / 3);
    Containers::arrayResize(phasefield, vertices.size());
}



void Viewer::drawShaderOptions(){
    if(!drawable) return;
    if(ImGui::TreeNode("Shader Options")) {
        ImGui::Text("Rendering: %3.2f FPS", Double(ImGui::GetIO().Framerate));
        static std::map<DrawableType, std::string> drawablemap{
                {DrawableType::MeshVisualizer, "Mesh Visualizer"},
                {DrawableType::PhongDiffuse,   "Phong Diffuse"},
                {DrawableType::FaceColored,   "Face Colored"},
                {DrawableType::FlatTextured,   "Flat Textured"}
        };

        if (ImGui::BeginCombo("##combodrawble", drawablemap[drawableType].c_str())) {
            for (auto const&[drtype, name] : drawablemap) {
                bool isSelected = (drawableType == drtype);
                if (ImGui::Selectable(name.c_str(), isSelected)) {
                    drawableType = drtype;
                    makeDrawableCurrent(drawableType);
                }
                if (isSelected)
                    ImGui::SetItemDefaultFocus();
            }
            ImGui::EndCombo();
        }

        if (drawableType == DrawableType::MeshVisualizer) {
            constexpr Float min = 0.f, max = 10.f;
            auto &d = dynamic_cast<MeshVisualizerDrawable &>(*drawable);
            ImGui::DragScalar("Wireframe Width", ImGuiDataType_Float, &d.wireframeWidth, .01f, &min, &max, "%f", 1);
            ImGui::ColorEdit3("Wireframe Color", d.wireframeColor.data());
            ImGui::ColorEdit3("Color", d.color.data());
            ImGui::DragScalar("Smoothness", ImGuiDataType_Float, &d.smoothness, .01f, &min, &max, "%f", 1);
        }

        static std::map<ColorMapType, std::string> colormapmap{
                {ColorMapType::Turbo,   "Turbo"},
                {ColorMapType::Magma,   "Magma"},
                {ColorMapType::Plasma,  "Plasma"},
                {ColorMapType::Inferno, "Inferno"},
                {ColorMapType::Viridis, "Viridis"}
        };

        if (drawableType == DrawableType::PhongDiffuse || drawableType == DrawableType::FlatTextured) {
            if (ImGui::BeginCombo("##combocmm", colormapmap[colorMapType].c_str())) {
                for (auto const&[cmtype, name] : colormapmap) {
                    bool isSelected = (colorMapType == cmtype);
                    if (ImGui::Selectable(name.c_str(), isSelected)) {
                        colorMapType = cmtype;
                        makeDrawableCurrent(drawableType);
                    }
                    if (isSelected)
                        ImGui::SetItemDefaultFocus();
                }
                ImGui::EndCombo();
            }
        }

        ImGui::TreePop();
    }
}

void Viewer::viewportEvent(ViewportEvent& event) {
    GL::defaultFramebuffer.setViewport({{}, event.framebufferSize()});
    if(camera)
        camera->reshape(event.windowSize(), event.framebufferSize());

    imgui.relayout(Vector2{event.windowSize()}/event.dpiScaling(),
                     event.windowSize(), event.framebufferSize());
}


void Viewer::keyPressEvent(KeyEvent& event) {
    if(imgui.handleKeyPressEvent(event)) {
        event.setAccepted();
        return;
    }

    if(!camera) return;

    switch(event.key()) {
        case KeyEvent::Key::L:
            if(camera->lagging() > 0.0f) {
                Debug{} << "Lagging disabled";
                camera->setLagging(0.0f);
            } else {
                Debug{} << "Lagging enabled";
                camera->setLagging(0.85f);
            }
            break;
        case KeyEvent::Key::R:
            camera->reset();
            break;
        case KeyEvent::Key::LeftCtrl :
            brushing = true;
            GL::defaultFramebuffer.mapForRead(GL::DefaultFramebuffer::ReadAttachment::Front);
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

    if(event.key() == Viewer::KeyEvent::Key::LeftCtrl){
        brushing = false;
        GL::defaultFramebuffer.mapForRead(GL::DefaultFramebuffer::ReadAttachment::None);
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

void Viewer::mousePressEvent(MouseEvent& event) {
    if(imgui.handleMousePressEvent(event)){
        event.setAccepted();
        return;
    }

    if(brushing) {
        point = Vector3d(unproject(event.position()));
        geodesicSearch();
        targetDist = 0;
        stopPainting = false;
        event.setAccepted();
        return;
    }

    if(event.button() == MouseEvent::Button::Middle){
            trackingMouse = true;
            ///* Enable mouse capture so the mouse can drag outside of the window */
            ///** @todo replace once https://github.com/mosra/magnum/pull/419 is in */
            SDL_CaptureMouse(SDL_TRUE);

            camera->initTransformation(event.position());

            event.setAccepted();
            redraw(); /* camera has changed, redraw! */
    }
}

void Viewer::mouseReleaseEvent(MouseEvent& event) {
    if(imgui.handleMouseReleaseEvent(event)) {
        event.setAccepted();
        return;
    }

    if(!stopPainting) {
        stopPainting = true;
        event.setAccepted();
        return;
    }

    if(event.button() == MouseEvent::Button::Middle) {
        /* Disable mouse capture again */
        /** @todo replace once https://github.com/mosra/magnum/pull/419 is in */
        if (trackingMouse) {
            SDL_CaptureMouse(SDL_FALSE);
            trackingMouse = false;
            event.setAccepted();
        }
    }
}

void Viewer::mouseMoveEvent(MouseMoveEvent& event) {
    if(imgui.handleMouseMoveEvent(event)){
        event.setAccepted();
        return;
    }

    //if(brushing){
    //    stopPainting = true;
    //    event.setAccepted();
    //    return;
    //}

    if(trackingMouse) {
        if (event.modifiers() & MouseMoveEvent::Modifier::Shift)
            camera->translate(event.position());
        else camera->rotate(event.position());

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

    if(!camera) return;

    const Float delta = event.offset().y();
    if(Math::abs(delta) < 1.0e-2f) return;

    camera->zoom(delta);

    event.setAccepted();
    redraw(); /* camera has changed, redraw! */
}


void Viewer::tickEvent()
{
    if(!stopPainting) {
        paint();
        reuploadVertices(vertexBuffer, meshData);
        redraw();
        return;
    }

    //exclusive events
    {
        std::lock_guard l(mutex);
        if (update & VisualizationFlag::ConnectedComponents) {
            faceTexture->setSubImage(0, {}, ImageView2D{Magnum::PixelFormat::RGB8Srgb, {(Int)faceColors.size(), 1}, faceColors});
            update &= ~VisualizationFlag::ConnectedComponents;
        }
        if (update & VisualizationFlag::GeodesicWeights) {
            faceTexture->setSubImage(0, {}, ImageView2D{Magnum::PixelFormat::RGB8Srgb, {(Int)faceColors.size(), 1}, faceColors});
            update &= ~VisualizationFlag::GeodesicWeights;
        }
        if (update & VisualizationFlag::Gradient) {
            reuploadVertices(vertexBuffer, meshData);
            update &= ~VisualizationFlag::Gradient;
        }
        if (update & VisualizationFlag::Phasefield) {
            reuploadVertices(vertexBuffer, meshData);
            update &= ~VisualizationFlag::Phasefield;
        }
        for (auto const& f : problem.functionals) {
            f->metaData->updateVis();
        }
    }

}

void Viewer::drawEvent() {
    GL::defaultFramebuffer.clear(GL::FramebufferClear::Color|GL::FramebufferClear::Depth);
    imgui.newFrame();

    /* Enable text input, if needed */
    if(ImGui::GetIO().WantTextInput && !isTextInputActive())
        startTextInput();
    else if(!ImGui::GetIO().WantTextInput && isTextInputActive())
        stopTextInput();

    /* draw scene */
    if(camera){
        bool camChanged = camera->update();
        camera->draw(drawableGroup);
    }

    //draw modifiers
    drawSubdivisionOptions();
    drawPrimitiveOptions();
    drawBrushOptions();
    drawOptimizationContext();
    drawShaderOptions();

    imgui.updateApplicationCursor(*this);

    /* Render ImGui window */
    {
        GL::Renderer::enable(GL::Renderer::Feature::Blending);
        GL::Renderer::disable(GL::Renderer::Feature::FaceCulling);
        GL::Renderer::disable(GL::Renderer::Feature::DepthTest);
        GL::Renderer::enable(GL::Renderer::Feature::ScissorTest);

        imgui.drawFrame();

        GL::Renderer::disable(GL::Renderer::Feature::ScissorTest);
        GL::Renderer::enable(GL::Renderer::Feature::DepthTest);
        GL::Renderer::enable(GL::Renderer::Feature::FaceCulling);
        GL::Renderer::disable(GL::Renderer::Feature::Blending);
    }

    swapBuffers();
    redraw();
}




