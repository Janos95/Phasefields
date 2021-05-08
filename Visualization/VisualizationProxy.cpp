//
// Created by janos on 09.06.20.
//

#include "VisualizationProxy.h"
#include "Viewer.h"
#include "Enums.h"
#include "C1Functions.h"

#include <Magnum/ImageView.h>
#include <Magnum/GL/TextureFormat.h>
#include <Magnum/GL/Sampler.h>
#include <Magnum/DebugTools/ColorMap.h>
#include <Magnum/Math/FunctionsBatch.h>

#include <Corrade/Utility/Algorithms.h>
#include <Corrade/Containers/GrowableArray.h>

#include <random>

namespace Phasefield {

using namespace Mg::Math::Literals;

const Color4 kelly_colors[] = {
    0xFFFFB300_rgbf, //Vivid Yellow
    0xFF803E75_rgbf, //Strong Purple
    0xFFFF6800_rgbf, //Vivid Orange
    0xFFA6BDD7_rgbf, //Very Light Blue
    0xFFC10020_rgbf, //Vivid Red
    0xFFCEA262_rgbf, //Grayish Yellow
    0xFF817066_rgbf, //Medium Gray
    0xFF007D34_rgbf, //Vivid Green
    0xFFF6768E_rgbf, //Strong Purplish Pink
    0xFF00538A_rgbf, //Strong Blue
    0xFFFF7A5C_rgbf, //Strong Yellowish Pink
    0xFF53377A_rgbf, //Strong Violet
    0xFFFF8E00_rgbf, //Vivid Orange Yellow
    0xFFB32851_rgbf, //Strong Purplish Red
    0xFFF4C800_rgbf, //Vivid Greenish Yellow
    0xFF7F180D_rgbf, //Strong Reddish Brown
    0xFF93AA00_rgbf, //Vivid Yellowish Green
    0xFF593315_rgbf, //Deep Yellowish Brown
    0xFFF13A13_rgbf, //Vivid Reddish Orange
    0xFF232C16_rgbf //Dark Olive Green
};

static Array<Color4> g_colors;

#define USE_KELLY 0

Array<Color4>& getColors(size_t n) {
    if(n == Invalid) {
        return g_colors;
    }
    if(n != g_colors.size()) {
        arrayResize(g_colors, n);
        if(n > 20 || !USE_KELLY) {
            Deg hue = 100.0_degf;
            Deg step = 360._degf/n;
            for(size_t i = 0; i < n; ++i) {
                g_colors[i] = Color4::fromHsv({hue += step, 0.9f, 0.7f});
            }
        } else {
            for(size_t i = 0; i < n; ++i) {
                g_colors[i] = kelly_colors[i];
            }
        }
    }
    return g_colors;
}

void optimizeColors(Array<size_t> const& neighbors, Array<size_t> const& starts) {
    size_t segmentCount = starts.size() - 1;
    auto& colors = getColors(segmentCount);

    std::random_device device;
    std::mt19937 gen(device());
    std::uniform_int_distribution<size_t> idxDist(0, segmentCount - 1);

    for(size_t i = 0; i < 1'000'000; ++i) {

        size_t indices[2] = {idxDist(gen), idxDist(gen)};
        float costOrig = 0;
        float costSwapped = 0;

        for(size_t l = 0; l < 2; ++l) {
            size_t idx = indices[l];
            Color4& colorOrig = colors[indices[l]];
            Color4& colorSwapped = colors[indices[l^1]];
            for(size_t j = starts[idx]; j < starts[idx+1]; ++j) {
                size_t k = neighbors[j];
                CORRADE_ASSERT(k < colors.size(), "Color Opt: Index out of bounds",);
                costOrig += (colorOrig - colors[k]).dot();
                costSwapped += (colorSwapped - colors[k]).dot();
            }
        }

        if(costSwapped > costOrig) {
            std::swap(colors[indices[0]], colors[indices[1]]);
        }
    }
}

VisualizationProxy::VisualizationProxy(Viewer& v) : viewer(v) {
    setDefaultCallback();
}

//void VisualizationProxy::setFaceColors(Containers::ArrayView<double>& data){
//    std::lock_guard l(mutex);
//}
//
//void setFaceColors(Containers::ArrayView<Mg::Color3ub>& data){
//
//}


void VisualizationProxy::upload() {
    if(updateVertexBuffer) {
        viewer.mesh.uploadVertexBuffer(viewer.vertexBuffer);
        updateVertexBuffer = false;
    }
}

void VisualizationProxy::drawSegmentation() {

    Tree& tree = viewer.tree;

    auto& colors = getColors(tree.nodeCountOnLevel(level)*2);
    auto vertexColors = viewer.mesh.colors();

    for(Color4& c : vertexColors) c = Color4{};

    size_t n = tree.vertexCount();
    SmootherStep smoothStep;

    tree.computeWeightsOfAncestorsOfLevel(level);

    size_t idx = 0;
    for(Node node : tree.nodesOnLevel(level)) {
        for(size_t i = 0; i < n; ++i) {
            vertexColors[i].rgb() += node.temporary()[i]*smoothStep.eval(node.phasefield()[i])*colors[2*idx].rgb();
            vertexColors[i].rgb() += node.temporary()[i]*smoothStep.eval(-node.phasefield()[i])*colors[2*idx + 1].rgb();
        }
        ++idx;
    }
    shaderConfig = ShaderConfig::VertexColors;
}

void VisualizationProxy::redraw() {
    drawCb(viewer.currentNode);
    viewer.redraw();
    updateVertexBuffer = true;
}

void VisualizationProxy::setDefaultCallback() {
    switch(option) {
        case VisOption::Segmentation :
            drawCb = [this](Node) { drawSegmentation(); };
            break;
        case VisOption::Phasefield :
            scale = 0.5;
            offset = 1;
            drawCb = [this](Node node) { drawValues(node.phasefield()); };
            break;
        case VisOption::Weight :
            scale = 1;
            offset = 0;
            drawCb = [this](Node node) { drawValues(node.temporary()); };
            break;
    }
    releaseCb = []{};
}

void VisualizationProxy::setCallbacks(UniqueFunction<void(Node)> draw, UniqueFunction<void()> release) {
    releaseCb();
    drawCb = std::move(draw);
    releaseCb = std::move(release);
}

void VisualizationProxy::drawValues(VertexDataView<double> const& values) {
    CORRADE_INTERNAL_ASSERT(values.size() == viewer.mesh.vertexCount());
    for(Vertex v : viewer.mesh.vertices()) {
        viewer.mesh.scalar(v) = scale*(values[v] + offset);
    }
    shaderConfig = ShaderConfig::ColorMaps;
}

void VisualizationProxy::drawValuesNormalized(VertexDataView<double> const& values) {
    if(!customMapping) {
        double min, max, w;
        std::tie(min, max) = Math::minmax(ArrayView<double>(values));
        if(max - min < 1e-10)
            w = 1;
        else
            w = max - min;
        Debug{} << "Min" << min << ", Max" << max;
        scale = 1/w;
        offset = -min;
    }
    drawValues(values);
}

}