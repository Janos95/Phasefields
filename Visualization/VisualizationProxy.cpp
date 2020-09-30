//
// Created by janos on 09.06.20.
//

#include "VisualizationProxy.h"
#include "Viewer.h"
#include "Enums.h"
#include "normalizeInto.hpp"
#include "C1Functions.h"

#include <Magnum/ImageView.h>
#include <Magnum/PixelFormat.h>
#include <Magnum/GL/TextureFormat.h>
#include <Magnum/GL/Sampler.h>
#include <Magnum/GL/Buffer.h>
#include <Magnum/GL/Mesh.h>
#include <Magnum/DebugTools/ColorMap.h>
#include <Magnum/Trade/MeshData.h>

#include <Corrade/Utility/Algorithms.h>
#include <Corrade/Containers/StaticArray.h>
#include <Corrade/Containers/GrowableArray.h>

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

Array<Color4>& getColors(size_t n) {
    static Array<Color4> colors;
    if(n != colors.size()) {
        arrayResize(colors, n);
        Deg hue = 100.0_degf;
        for(auto& c : colors) /* generate random colors */
            c = Color4::fromHsv({hue += 137.5_degf, 0.9f, 0.9f});
    }
    return colors;
}


VisualizationProxy::VisualizationProxy(Viewer& v) : viewer(v) {
    using L = std::initializer_list<std::pair<const char*, StaticArrayView<256, const Vector3ub>>>;
    for(auto [name, colors] : L{
            {"Turbo",   Magnum::DebugTools::ColorMap::turbo()},
            {"Magma",   Magnum::DebugTools::ColorMap::magma()},
            {"Plasma",  Magnum::DebugTools::ColorMap::plasma()},
            {"Inferno", Magnum::DebugTools::ColorMap::inferno()},
            {"Viridis", Magnum::DebugTools::ColorMap::viridis()}
    }) {
        Array<Color4> converted{256};
        for(size_t i = 0; i < 256; ++i) {
            converted[i] = Color4{Math::unpack<Color3>(colors[i])};
        }
        arrayAppend(colorMaps, InPlaceInit, name, std::move(converted));
    }
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

    auto& colors = getColors(tree.numLeafs*2);
    auto vertexColors = viewer.mesh.colors();

    for(double& w : tree.root().temporary()) w = 1.;
    for(Color4& c : vertexColors) c = Color4{};

    size_t n = tree.vertexCount();
    SmootherStep smoothStep;

    //auto [min,max] = Math::minmax(phasefields[0]);
    //Debug{} << min << " " << max;

    for(Node node : tree.internalNodes()) {
        Debug{} << node;
        auto weights = node.temporary();
        for(size_t i = 0; i < n; ++i) {
            if(node.hasLeftChild()) {
                Node leftChild = node.leftChild();
                leftChild.temporary()[i] = smoothStep.eval(node.phasefield()[i])*weights[i];
            }
            if(node.hasRightChild()) {
                Node rightChild = node.rightChild();
                rightChild.temporary()[i] = smoothStep.eval(-node.phasefield()[i])*weights[i];
            }
        }
    }

    size_t leafIdx = 0;
    for(Node leaf : tree.leafs()) {
        for(size_t i = 0; i < n; ++i) {
            vertexColors[i].rgb() += leaf.temporary()[i]*smoothStep.eval(leaf.phasefield()[i])*colors[2*leafIdx].rgb();
            vertexColors[i].rgb() += leaf.temporary()[i]*smoothStep.eval(-leaf.phasefield()[i])*colors[2*leafIdx + 1].rgb();
        }
        ++leafIdx;
    }

}

void VisualizationProxy::redraw() {
    drawCb(&viewer);
    viewer.redraw();
    updateVertexBuffer = true;
}

void VisualizationProxy::setDefaultCallback() {
    isDefaultCallback = true;
    switch(option) {
        case VisOption::Segmentation :
            drawCb = [this](Viewer*) { drawSegmentation(); };
            break;
        case VisOption::Phasefield :
            drawCb = [this](Viewer*) { drawValues(viewer.currentNode.phasefield(), [](double x){ return 0.5*(x + 1); }); };
            break;
        case VisOption::Weight :
            drawCb = [this](Viewer*) { drawValues(viewer.currentNode.temporary(), [](double x){ return x; }); };
            break;
    }
    releaseCb = []{};
}

void VisualizationProxy::setCallbacks(UniqueFunction<void(Viewer*)> draw, UniqueFunction<void()> release) {
    releaseCb();
    drawCb = std::move(draw);
    releaseCb = std::move(release);
    isDefaultCallback = false;
}

void VisualizationProxy::drawValues(VertexDataView<double> values, UniqueFunction<double(double)> tf) {
    auto& map = colorMaps[0].second;
    CORRADE_INTERNAL_ASSERT(values.size() == viewer.mesh.vertexCount());
    for(Vertex v : viewer.mesh.vertices()) {
        double i = tf(values[v])*255.;
        Color4 c;
        if(i <= 0.) {
            viewer.mesh.color(v) = map[0];
        }
        else if(i >= 255.) {
            viewer.mesh.color(v) = map[255];
        }
        else {
            double u = Math::ceil(i);
            double l = Math::floor(i);
            viewer.mesh.color(v) = (u - i)*map[size_t(l)] + (i - l)*map[size_t(u)];
        }
    }
}

}