//
// Created by janos on 13.03.20.
//

#include "upload.hpp"

#include <Magnum/GL/TextureFormat.h>
#include <Magnum/GL/Sampler.h>
#include <Magnum/MeshTools/Compile.h>

using namespace Magnum;
using namespace Corrade;

Object upload(Trade::MeshData& meshData, ImageView2D image, CompileFlag flags) {
    Object object;
    object.textureDiffuse.reset(new GL::Texture2D{});
    object.textureDiffuse->setMagnificationFilter(GL::SamplerFilter::Linear)
            .setMinificationFilter(GL::SamplerFilter::Linear, GL::SamplerMipmap::Linear)
            .setWrapping(GL::SamplerWrapping::ClampToEdge)
            .setMaxAnisotropy(GL::Sampler::maxMaxAnisotropy())
            .setStorage(Math::log2(4096) + 1, GL::TextureFormat::RGBA8, {4096, 4096})
            .setSubImage(0, {}, image)
            .generateMipmap();
    GL::Buffer indices, vertices;
    indices.setData(meshData.indexData());
    vertices.setData(meshData.vertexData());
    object.vertices = std::move(vertices);
    object.indices = std::move(indices);
    object.mesh = MeshTools::compile(meshData, MeshTools::CompileFlag::GenerateFlatNormals);

    return object;
}
