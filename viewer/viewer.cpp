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
#include <Magnum/Math/Matrix4.h>
#include <Magnum/GL/PixelFormat.h>
#include <Magnum/Math/FunctionsBatch.h>
#include <Magnum/ImGuiIntegration/Context.hpp>
#include <Magnum/Trade/MeshData.h>
#include <Magnum/Shaders/Flat.h>
#include <Magnum/Image.h>
#include <Magnum/GL/TextureFormat.h>
#include <Magnum/PixelFormat.h>
#include <Magnum/ImageView.h>
#include <Magnum/GL/Renderer.h>
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




Viewer::Viewer(int argc, char** argv):
    Platform::Application{{argc,argv},NoCreate},
    dirichletScaling(1.), doubleWellScaling(1.), connectednessScaling(1.),
    phasefieldCallback{*this}
{
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



    paths = new Paths(&scene, drawableGroup);
    object = new Object3D(&scene);

    GL::Renderer::enable(GL::Renderer::Feature::DepthTest);
    GL::Renderer::enable(GL::Renderer::Feature::FaceCulling);

    /* Start the timer, loop at 60 Hz max */
    setSwapInterval(1);
    setMinimalLoopPeriod(16);
}

solver::Status UpdatePhasefield::operator ()(solver::IterationSummary const&){
    std::lock_guard l(viewer.mutex);
    if(viewer.visFlags & VisualizationFlag::Phasefield)
        Utility::copy(viewer.parameters, viewer.phasefield);
    viewer.update |= VisualizationFlag::Phasefield;
    return solver::Status::CONTINUE;
}

void Viewer::drawSubdivisionOptions() {
    if (ImGui::TreeNode("Subdivisions")) {
        ImGui::Text("Currently we have %d vertices and %d faces", (int)vertices.size(), (int)indices.size() / 3);
        constexpr int step = 1;
        static std::uint32_t subs = 0;
        ImGui::InputScalar("Number of Sibdivision (wrt. orignal mesh)", ImGuiDataType_U32, &subs, &step, nullptr, "%d");
        if (ImGui::Button("do subdivision")) {
            if(subs >= numSubdivisions) subdivide(subs - numSubdivisions, meshData, phasefield, vertices, indices, meshData);
            else subdivide(subs, original, phasefield, vertices, indices, meshData);
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
            d->texture = visFlags & VisualizationFlag::ConnectedComponents ? componentsTexture.get() : wsTexture.get();
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
        updateIntenalDataStrucutres();
        upload(mesh, vertexBuffer, indexBuffer, meshData);
        updateFunctionals(problem.functionals);
        updateFunctionals(problem.constraints);

        makeDrawableCurrent(drawableType);
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


void Viewer::dumpToPly(std::string const& path, ConnectednessMetaData<Double> const& meta){
    auto triangles = Cr::Containers::arrayCast<Vector3ui>(indices);
    std::ofstream out("/tmp/phase.ply");
    out << "ply" << std::endl;
    out << "format ascii 1.0\n";
    out << "element vertex " << vertices.size() << '\n';
    out << "property float x\n";
    out << "property float y\n";
    out << "property float z\n";
    out << "property float j\n";
    out << "property float u\n";
    out << "element face " << triangles.size() << '\n';
    out << "property list uchar int vertex_indices\n";
    out << "property int c\n";
    out << "property float w\n";
    out << "end_header\n";

    for (int i = 0; i < vertices.size(); ++i) {
        for (int j = 0; j < 3; ++j) {
            out << vertices[i][j] << ' ';
        }
        out << gradient[i] << ' ' << phasefield[i] << '\n';
    }

    for (int i = 0; i < triangles.size(); ++i) {
        out << "3 ";
        for (int j = 0; j < 3; ++j) {
            out << triangles[i][j] << ' ';
        }
        out << components[i] << ' ' << ws[i] << '\n';
    }
}

/**
 * returns true if functional is now handling exclusive visualizations
 * Also note that we can safely read from meta.flags since the optimization thread
 * only reads from those as well (even taking the lock in case we modify)
 */
bool Viewer::drawConnectednessConstraintOptions(ConnectednessConstraint<Double>& f, VisualizationFlags& fl){
    auto& meta = f.getMetaData();
    double a = meta.a, b = meta.b;
    bool makeExclusive = false;
    if(dragDoubleRange2("Positive Interval", &a, &b, 0.01f, -1.f, 1.f, "Min: %.2f", "Max: %.2f",1.f)){
        std::lock_guard l(mutex);
        meta.a = a; meta.b = b;
    }
    bool visPaths = static_cast<bool>(meta.flags & VisualizationFlag::Paths);
    if(ImGui::Checkbox("Visualize Shortest Paths", &visPaths)){
        std::lock_guard l(mutex);
        if(visPaths) {
            paths->muted = false;
            meta.flags |= VisualizationFlag::Paths;
        }
        else {
            meta.flags &= ~VisualizationFlag::Paths;
            paths->muted = true;
        }
    }
    ImGui::SameLine();
    bool visComponents = static_cast<bool>(meta.flags & VisualizationFlag::ConnectedComponents);
    if(ImGui::Checkbox("Connected Components", &visComponents)){
        std::lock_guard l(mutex);
        if(visComponents) {
            makeDrawableCurrent(DrawableType::FaceColored);
            meta.flags &= NonExclusiveFlags; // deselect all
            meta.flags |= VisualizationFlag::ConnectedComponents;
            update |= VisualizationFlag::ConnectedComponents;
            makeExclusive = true;
        }
        else meta.flags &= ~VisualizationFlag::ConnectedComponents;
    }
    bool visGeodesicWeights = static_cast<bool>(meta.flags & VisualizationFlag::GeodesicWeights);
    if(ImGui::Checkbox("Geodesic Weights", &visGeodesicWeights)){
        std::lock_guard l(mutex);
        if(visGeodesicWeights) {
            makeDrawableCurrent(DrawableType::FaceColored);
            meta.flags &= NonExclusiveFlags; // deselect all
            meta.flags |= VisualizationFlag::GeodesicWeights;
            update |= VisualizationFlag::GeodesicWeights;
            makeExclusive = true;
        }
        else meta.flags &= ~VisualizationFlag::GeodesicWeights;
    }

    ImGui::SameLine();
    bool visGradient = static_cast<bool>(meta.flags & VisualizationFlag::Gradient);
    if(ImGui::Checkbox("Gradient", &visGradient)){
        std::lock_guard l(mutex);
        if(visGradient) {
            makeDrawableCurrent(DrawableType::PhongDiffuse);
            meta.flags &= NonExclusiveFlags; // deselect all
            meta.flags |= VisualizationFlag::Gradient;
            update |= VisualizationFlag::Gradient;
            makeExclusive = true;
        }
        else meta.flags &= ~VisualizationFlag::Gradient;
    }
    fl = meta.flags;
    return makeExclusive;
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
        auto& fs = problem.functionals;
        for(int i = 0; i < fs.size(); ++i){
            auto& f = *fs[i];
            ImGui::PushID(++nodeCount);
            ImGui::Separator();

            if(ImGui::Button("Remove")) toRemove = i;
            ImGui::SameLine();

            switch(f.type){
                case FunctionalType::DoubleWellPotential :
                    ImGui::Text("Double Well Potential");
                    break;
                case FunctionalType::DirichletEnergy :
                    ImGui::Text("Dirichlet Energy");
                    break;
                case FunctionalType::Area1 :
                    ImGui::Text("Area Regularization Quadratic");
                    drawAreaRegularization1Options(dynamic_cast<AreaRegularizer1&>(f));
                    break;
                case FunctionalType::Area2 :
                    ImGui::Text("Area Regularization Smooth Step");
                    drawAreaRegularization2Options(dynamic_cast<AreaRegularizer2&>(f));
                    break;
                case FunctionalType::Connectedness :
                    ImGui::Text("Connectedness Constraint");
                    VisualizationFlags flags;
                    if(drawConnectednessConstraintOptions(dynamic_cast<ConnectednessConstraint<Double>&>(f), flags)){
                        auto funcPtr = fs[i].get();
                        if(exclusiveVisualizer && exclusiveVisualizer != funcPtr) {
                            std::lock_guard l(mutex);
                            exclusiveVisualizer->metaData->flags &= NonExclusiveFlags;
                        }
                        exclusiveVisualizer = fs[i].get();
                        visFlags &= flags | NonExclusiveFlags; //set exclusive flags from functional
                    } else visFlags |= flags;
                    break;
            }

            drawLoss(f.metaData->loss, nodeCount);
            ImGui::PopID();
        }

        if(nodeCount) ImGui::Separator();

        if(toRemove >= 0){
            auto f = problem.functionals[toRemove].get();
            visFlags &= ~f->metaData->flags; //clear everything it was displaying
            if(f == exclusiveVisualizer)
                exclusiveVisualizer = nullptr;

            std::swap(fs[toRemove], fs.back());
            Containers::arrayResize(fs, fs.size() - 1);
        }

        bool visPhasefield = static_cast<bool>(visFlags & VisualizationFlag::Phasefield);
        if(ImGui::Checkbox("Draw Phasefield", &visPhasefield)){
            if(visPhasefield){
                makeDrawableCurrent(DrawableType::PhongDiffuse);
                auto textureCoords = meshData.mutableAttribute<Vector2>(Trade::MeshAttribute::TextureCoordinates);
                for (int i = 0; i < phasefield.size(); ++i)
                    textureCoords[i].x() = .5f * (phasefield[i] + 1.f);
                visFlags &= NonExclusiveFlags;
                visFlags |= VisualizationFlag::Phasefield;
                update |= VisualizationFlag::Phasefield;
                if(exclusiveVisualizer){
                    std::lock_guard l(mutex);
                    exclusiveVisualizer->metaData->flags &= NonExclusiveFlags;
                    exclusiveVisualizer = nullptr;
                }
            } else visFlags &= ~VisualizationFlag::Phasefield;
        }

        static Double epsilon = 0.1;
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

    Containers::arrayResize(parameters, phasefield.size());
    Utility::copy(phasefield, parameters);

    Containers::Array<solver::iteration_callback_type> callbacks;
    for(auto& f : problem.functionals)
        Containers::arrayAppend(callbacks, Containers::InPlaceInit, *f->metaData);
    Containers::arrayAppend(callbacks, Containers::InPlaceInit, phasefieldCallback);
    options.callbacks = std::move(callbacks);

    optimizing = true;
    thread = std::thread([this]{
        solver::Summary summary;
        solve(options, problem, parameters.data(), &summary);
        //Debug{} << summary.briefReport().c_str();
    });
}

void Viewer::stopOptimization() {
    if(optimizing.exchange(false))
        thread.join();
}

Containers::Pointer<Functional> Viewer::makeFunctional(FunctionalType type) {
    auto ts = Containers::arrayCast<Vector3ui>(indices);
    switch (type) {
        case FunctionalType::Area1 : {
            auto p = Containers::pointer<AreaRegularizer1>(vertices, ts);
            p->type = FunctionalType::Area1; //@todo somehow this is not working
            return p;
        }
        case FunctionalType::Area2 : {
            auto p = Containers::pointer<AreaRegularizer2>(vertices, ts);
            p->type = FunctionalType::Area2;
            return p;
        }
        case FunctionalType::DirichletEnergy : {
            auto p = Containers::pointer<DirichletEnergy>(vertices, ts);
            p->metaData->scaling = dirichletScaling;
            return p;
        }
        case FunctionalType::DoubleWellPotential : {
            auto p = Containers::pointer<DoubleWellPotential>(vertices, ts);
            p->metaData->scaling = doubleWellScaling;
            p->type = FunctionalType::DoubleWellPotential;
            return p;
        }
        case FunctionalType::Connectedness :
            auto p = Containers::pointer<ConnectednessConstraint<Double>>(vertices, ts);
            p->metaData->scaling = connectednessScaling;
            p->getMetaData().viewer = this;
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

void Viewer::updateIntenalDataStrucutres(){
    meshData = preprocess(original, CompileFlag::GenerateSmoothNormals|CompileFlag::AddTextureCoordinates);
    CORRADE_ASSERT(meshData.hasAttribute(Trade::MeshAttribute::TextureCoordinates),"error", );
    auto points = meshData.attribute<Vector3>(Trade::MeshAttribute::Position);
    Containers::arrayResize(vertices, points.size());
    for (int i = 0; i < points.size(); ++i)
        vertices[i] = Vector3d(points[i]);
    indices = meshData.indicesAsArray();
    Containers::arrayResize(phasefield, vertices.size());
    upload(mesh, vertexBuffer, indexBuffer, meshData);
    const Magnum::Vector2i size(indices.size() / 3, 1);
    wsTexture.reset(new GL::Texture2D{});
    wsTexture->setStorage(1, Magnum::GL::TextureFormat::RGB8, size);
    componentsTexture.reset(new GL::Texture2D{});
    componentsTexture->setStorage(1, Magnum::GL::TextureFormat::RGB8, size);
    Containers::arrayResize(gradient, vertices.size());
    Containers::arrayResize(ws, indices.size() / 3);
    Containers::arrayResize(components, indices.size() / 3);
    numComponents = 0;
    update = visFlags; //update everything
}


void Viewer::updateFaceColorTextureWithComponents() {
    Deg hue = 42.0_degf;
    Containers::Array<Color3ub> randomColors(Containers::NoInit, numComponents);
    for (int i = 0; i < numComponents; ++i)
        randomColors[i] = Color3ub::fromHsv({hue += 137.5_degf, 0.75f, 0.9f});

    Containers::Array<Vector3ub> colors(Containers::NoInit, components.size());
    for (int i = 0; i < components.size(); ++i) {
        if (components[i] >= 0 && !randomColors.empty())
            colors[i] = randomColors[components[i]];
        else
            colors[i] = Color3ub(255,253,208);
    }

    const Magnum::Vector2i size(components.size(), 1);
    componentsTexture->setMinificationFilter(Magnum::SamplerFilter::Nearest)
                      .setMagnificationFilter(Magnum::SamplerFilter::Nearest)
                      .setWrapping(Magnum::SamplerWrapping::ClampToEdge) // or Repeat
                      .setSubImage(0, {}, ImageView2D{Magnum::PixelFormat::RGB8Srgb, size, colors});

}

void Viewer::updateFaceColorTextureWithWeights() {

    Containers::Array<Vector3ub> colors(Containers::NoInit, ws.size());
    auto colorMap = DebugTools::ColorMap::turbo();
    auto [min,max] = Math::minmax(ws);
    Double length = max - min;
    if(length < std::numeric_limits<Double>::epsilon())
        std::fill(colors.begin(), colors.end(), Color3ub(255,253,208));
    else{
        for (int i = 0; i < ws.size(); ++i) {
            int idx = static_cast<int>((ws[i] - min) / length * (colorMap.size() - 1.));
            CORRADE_ASSERT(0 <= idx && idx <= 255, "Bad index for colormap", );
            colors[i] = colorMap[idx];
        }
    }
    const Magnum::Vector2i size(ws.size(), 1);
    wsTexture->setMinificationFilter(Magnum::SamplerFilter::Nearest)
              .setMagnificationFilter(Magnum::SamplerFilter::Nearest)
              .setWrapping(Magnum::SamplerWrapping::ClampToEdge) // or Repeat
              .setSubImage(0, {}, ImageView2D{Magnum::PixelFormat::RGB8Srgb, size, colors});
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

    if(!camera) return;

    trackingMouse = true;
    ///* Enable mouse capture so the mouse can drag outside of the window */
    ///** @todo replace once https://github.com/mosra/magnum/pull/419 is in */
    SDL_CaptureMouse(SDL_TRUE);

    camera->initTransformation(event.position());

    event.setAccepted();
    redraw(); /* camera has changed, redraw! */
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

    if(!camera) return;
    /* Disable mouse capture again */
    /** @todo replace once https://github.com/mosra/magnum/pull/419 is in */
    if(trackingMouse){
        SDL_CaptureMouse(SDL_FALSE);
        trackingMouse = false;
        event.setAccepted();
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

    if(!camera || !event.buttons()) return;

    if(event.modifiers() & MouseMoveEvent::Modifier::Shift)
        camera->translate(event.position());
    else camera->rotate(event.position());

    event.setAccepted();
    redraw(); /* camera has changed, redraw! */
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

    {
        std::lock_guard l(mutex);
        if (update & visFlags & VisualizationFlag::Paths) {
            Containers::arrayResize(paths->instanceData, Containers::NoInit, instanceData.size());
            Utility::copy(instanceData, paths->instanceData);
            update &= ~VisualizationFlag::Paths;
        }
        if (update & visFlags & VisualizationFlag::ConnectedComponents) {
            auto& f = dynamic_cast<ConnectednessConstraint<Double>&>(*exclusiveVisualizer);
            updateFaceColorTextureWithComponents();
            auto& d = dynamic_cast<FaceColorDrawable&>(*drawable);
            update &= ~VisualizationFlag::ConnectedComponents;
        }
        if (update & visFlags & VisualizationFlag::GeodesicWeights) {
            auto& f = dynamic_cast<ConnectednessConstraint<Double>&>(*exclusiveVisualizer);
            updateFaceColorTextureWithWeights();
            update &= ~VisualizationFlag::GeodesicWeights;
        }
        if (update & visFlags & VisualizationFlag::Gradient) {
            auto textureCoords = meshData.mutableAttribute<Vector2>(Trade::MeshAttribute::TextureCoordinates);
            for (int i = 0; i < gradient.size(); ++i)
                textureCoords[i].x() = gradient[i]; /* @note gradient needs to be normalized */
            reuploadVertices(vertexBuffer, meshData);
            update &= ~VisualizationFlag::Gradient;
        }
        if (update & visFlags & VisualizationFlag::Phasefield) {
            auto textureCoords = meshData.mutableAttribute<Vector2>(Trade::MeshAttribute::TextureCoordinates);
            for (int i = 0; i < phasefield.size(); ++i)
                textureCoords[i].x() = .5f * (phasefield[i] + 1.f);
            reuploadVertices(vertexBuffer, meshData);
            update &= ~VisualizationFlag::Phasefield;
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




