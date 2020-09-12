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

using namespace Magnum;
using namespace Corrade;

VisualizationProxy::VisualizationProxy(Viewer& v) : viewer(v) {
}

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

void VisualizationProxy::setVertexColors(Containers::StridedArrayView1D<double> const& data) {
    CORRADE_INTERNAL_ASSERT(data.size() == viewer.mesh.vertexCount());
    for(size_t i = 0; i < viewer.mesh.vertexCount(); ++i)
        viewer.mesh.scalars()[i] = data[i];
    shaderConfig = ShaderConfig::ColorMaps;
    updateVertexBuffer = true;
}

void VisualizationProxy::setTag(Mg::Int tag) {
    activeTag = tag;
}

bool VisualizationProxy::isActiveTag(Mg::Int tag) const {
    return tag == activeTag;
}

void VisualizationProxy::setVertexColors(Tree& tree) {

    if(colors.size() != tree.numLeafs*2) {
        Containers::arrayResize(colors, Containers::NoInit, tree.numLeafs*2);
        Deg hue = 42.0_degf;
        for(auto& c : colors) /* genereate random colors */
            c = Color4::fromHsv({hue += 137.5_degf, 0.75f, 0.9f});
    }

    auto vertexColors = viewer.mesh.colors();
    for(auto& c: vertexColors) c = Color4{0};

    auto phasefields = tree.phasefields();
    auto prefixes = tree.temporaryData();
    size_t n = tree.vertexCount();
    SmootherStep smoothStep;

    //auto [min,max] = Math::minmax(phasefields[0]);
    //Debug{} << min << " " << max;

    size_t leafIdx = 0;
    tree.traverse([&](Node& node) -> bool {
        if(!node.isLeaf()) {
            for(size_t i = 0; i < n; ++i) {
                if(node.leftChild != Invalid) {
                    prefixes[node.leftChild][i] = smoothStep.eval(phasefields[node.idx][i]);
                }
                if(node.rightChild != Invalid) {
                    prefixes[node.rightChild][i] = smoothStep.eval(-phasefields[node.idx][i]);
                }
            }
        } else {
            for(size_t i = 0; i < n; ++i) {
                vertexColors[i].rgb() += prefixes[node.idx][i]*smoothStep.eval(phasefields[node.idx][i])*colors[2*leafIdx].rgb();
                vertexColors[i].rgb() += prefixes[node.idx][i]*smoothStep.eval(-phasefields[node.idx][i])*colors[2*leafIdx + 1].rgb();
            }
            ++leafIdx;
        }
        return true;
    });

    shaderConfig = ShaderConfig::VertexColors;
    updateVertexBuffer = true;
}

}