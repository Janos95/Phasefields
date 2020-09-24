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
#include <Magnum/Trade/MeshData.h>

#include <Corrade/Utility/Algorithms.h>
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

VisualizationProxy::VisualizationProxy(Viewer& v) : viewer(v) { setDefaultCallback(); }

//void VisualizationProxy::setFaceColors(Containers::ArrayView<double>& data){
//    std::lock_guard l(mutex);
//}
//
//void setFaceColors(Containers::ArrayView<Mg::Color3ub>& data){
//
//}


void VisualizationProxy::upload() {
    if(updateFaceTexture) {
        //viewer.faceTexture.setSubImage(0, {}, Mg::ImageView2D{Magnum::PixelFormat::RGB8Srgb, {(int)faceColors.size(), 1}, faceColors});
        //updateFaceTexture = false;
    }
    if(updateVertexBuffer) {
        viewer.mesh.uploadVertexBuffer(viewer.vertexBuffer);
        updateVertexBuffer = false;
    }
}

void VisualizationProxy::drawPhasefield() {
    auto data = viewer.currentNode.phasefield();
    CORRADE_INTERNAL_ASSERT(data.size() == viewer.mesh.vertexCount());
    for(size_t i = 0; i < viewer.mesh.vertexCount(); ++i)
        viewer.mesh.scalars()[i] = 0.5*(data[i] + 1.);
}

void VisualizationProxy::drawWeights() {
    auto data = viewer.currentNode.temporary();
    for(size_t i = 0; i < viewer.mesh.vertexCount(); ++i)
        viewer.mesh.scalars()[i] = 0.5*(data[i] + 1.);
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
    cb(&viewer);
    viewer.redraw();
    updateVertexBuffer = true;
}

void VisualizationProxy::setDefaultCallback() {
    isDefaultCallback = true;
    switch(option) {
        case VisOption::Segmentation :
            cb = [this](Viewer*) { drawSegmentation(); };
            shaderConfig = ShaderConfig::VertexColors;
            break;
        case VisOption::Phasefield :
            cb = [this](Viewer*) { drawPhasefield(); };
            shaderConfig = ShaderConfig::ColorMaps;
            break;
        case VisOption::Weight :
            cb = [this](Viewer*) { drawWeights(); };
            shaderConfig = ShaderConfig::ColorMaps;
            break;
    }
}

void VisualizationProxy::setCallback(UniqueFunction<void(Viewer*)>&& cb_) {
    cb = std::move(cb_);
    isDefaultCallback = false;
}


}