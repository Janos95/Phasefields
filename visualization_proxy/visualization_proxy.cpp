//
// Created by janos on 09.06.20.
//

#include "visualization_proxy.hpp"
#include "viewer.hpp"
#include "types.hpp"

#include <Magnum/ImageView.h>
#include <Magnum/PixelFormat.h>
#include <Magnum/GL/TextureFormat.h>
#include <Magnum/GL/Sampler.h>
#include <Magnum/GL/Buffer.h>
#include <Magnum/GL/Mesh.h>
#include <Magnum/Trade/MeshData.h>

#include <mutex>

using namespace Magnum;
using namespace Corrade;

std::mutex mutex;

VisualizationProxy::VisualizationProxy(Viewer& v) : viewer(v)
{
}

void VisualizationProxy::setFaceColors(Containers::ArrayView<double>& data){
    std::lock_guard l(mutex);
}

void setFaceColors(Containers::ArrayView<Mg::Color3ub>& data){

}
void setVertexColors(Containers::ArrayView<double>& data){

}

void VisualizationProxy::resize(std::size_t n){
    Containers::arrayResize(faceColors, Containers::NoInit, n);
}

void VisualizationProxy::update(){
    std::lock_guard l(mutex);
    if(updateFaceTexture){
        viewer.faceTexture.setSubImage(0, {}, Mg::ImageView2D{Magnum::PixelFormat::RGB8Srgb, {(int)faceColors.size(), 1}, faceColors});
        updateFaceTexture = false;
    }
    if(updateVertexUvs){
        auto data = viewer.vertexBuffer.map(0,
                                            viewer.vertexBuffer.size(),
                                            GL::Buffer::MapFlag::Write);

        CORRADE_CONSTEXPR_ASSERT(data, "could not map vertex data");
        Utility::copy(viewer.meshData.vertexData(), data);
        viewer.vertexBuffer.unmap();
        updateVertexUvs = false;
    }
}
