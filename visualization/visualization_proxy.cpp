//
// Created by janos on 09.06.20.
//

#include "visualization_proxy.hpp"
#include "viewer.hpp"
#include "types.hpp"
#include "normalizeInto.hpp"

#include <Magnum/ImageView.h>
#include <Magnum/PixelFormat.h>
#include <Magnum/GL/TextureFormat.h>
#include <Magnum/GL/Sampler.h>
#include <Magnum/GL/Buffer.h>
#include <Magnum/GL/Mesh.h>
#include <Magnum/Trade/MeshData.h>

#include <Corrade/Utility/Algorithms.h>
#include <Corrade/Containers/GrowableArray.h>

#include <mutex>

using namespace Magnum;
using namespace Corrade;

std::mutex mutex;

VisualizationProxy::VisualizationProxy(Viewer& v) : viewer(v)
{
}

//void VisualizationProxy::setFaceColors(Containers::ArrayView<double>& data){
//    std::lock_guard l(mutex);
//}
//
//void setFaceColors(Containers::ArrayView<Mg::Color3ub>& data){
//
//}


void VisualizationProxy::update(){
    std::lock_guard l(mutex);
    if(updateFaceTexture){
        //viewer.faceTexture.setSubImage(0, {}, Mg::ImageView2D{Magnum::PixelFormat::RGB8Srgb, {(int)faceColors.size(), 1}, faceColors});
        //updateFaceTexture = false;
    }
    if(updateVertexBuffer){
        auto data = viewer.vertexBuffer.map(0,
                                            viewer.vertexBuffer.size(),
                                            GL::Buffer::MapFlag::Write);

        CORRADE_CONSTEXPR_ASSERT(data, "could not map vertex data");
        Utility::copy(viewer.meshData.vertexData(), data);
        viewer.vertexBuffer.unmap();
        updateVertexBuffer = false;
    }

}

void VisualizationProxy::setVertexColors(Containers::ArrayView<double>& data){
    std::lock_guard l(mutex);
    auto coords = viewer.meshData.mutableAttribute<Vector2>(Trade::MeshAttribute::TextureCoordinates);
    normalizeInto(data, Containers::arrayCast<2, float>(coords).slice<1>());
    updateVertexBuffer = true;
}

void VisualizationProxy::setVertexColors(PhasefieldTree& tree) {
    std::lock_guard l(mutex);

    if(colors.size() != tree.numLeafs*2){
        Containers::arrayResize(colors, Containers::NoInit, tree.numLeafs*2);
        Deg hue = 42.0_degf;
        for(auto& c : colors)
            c = Color4::fromHsv({hue += 137.5_degf, 0.75f, 0.9f});
    }

    auto vertexColors = viewer.meshData.mutableAttribute<Color4>(Trade::MeshAttribute::Color);
    for(auto& c: vertexColors) c = Color4{0};

    auto phasefields = tree.phasefields();
    auto prefixes = tree.temps();
    auto n = tree.phasefieldSize;
    SmootherStep smoothStep;

    int leafIdx = 0;
    tree.traverse([&](PhasefieldNode& node) -> void {
        if(!node.isLeaf()){
            for(UnsignedInt i = 0; i < n; ++i){
                if(node.leftChild != PhasefieldNode::None){
                    double pos = smoothStep.eval(phasefields[node.leftChild][i]);
                    prefixes[node.leftChild][i] = pos*phasefields[node.idx][i];
                }
                if(node.rightChild != PhasefieldNode::None){
                    double neg = smoothStep.eval(-phasefields[node.leftChild][i]);
                    prefixes[node.rightChild][i] = neg*phasefields[node.idx][i];
                }
            }
        } else {
            for(UnsignedInt i = 0; i < n; ++i){
                vertexColors[i].rgb() += prefixes[node.idx][i] * colors[leafIdx].rgb();
            }
            ++leafIdx;
        }
    });

    updateVertexBuffer = true;
}
