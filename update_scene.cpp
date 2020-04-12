//
// Created by janos on 03.04.20.
//

#include "update_scene.hpp"
#include "upload.hpp"

#include <Corrade/Utility/Algorithms.h>

void UpdateScene::tickEvent(Scene& scene)  {
    std::lock_guard l(phasefieldData.mutex);
    switch (phasefieldData.status) {
        case PhasefieldData::Status::NothingChanged : return;
        case PhasefieldData::Status::PhasefieldUpdated: reuploadVertices(); break;
        case PhasefieldData::Status::Subdivided:
        case PhasefieldData::Status::NewMesh: makeScene(scene); break;
    }
    phasefieldData.status = PhasefieldData::Status::NothingChanged;
    scene.setDirty();
}

void UpdateScene::makeScene(Scene& scene){
    //scene.reset(); //@todo note sure if the reset function correctly cleans up the drawable group

    //for(auto& v : md.mutableAttribute<Vector2>(md.attributeId(Trade::MeshAttribute::TextureCoordinates)))
    //        v = Vector2{1.f};
    //for(auto v : phasefieldData.meshData.textureCoordinates2DAsArray())
    //    Debug{} << v;

    //Debug{} << "is indexed " << phasefieldData.meshData.isIndexed();

    upload(phasefieldData, phasefieldData.meshData);
    phasefieldData.drawable = scene.addNode("mesh", phasefieldData, DrawableType::ColorMapPhong);
    phasefieldData.type = DrawableType::ColorMapPhong;
}

void UpdateScene::reuploadVertices()
{
    GL::Buffer& vertices = phasefieldData.vertices;
    auto data = vertices.map(0,
                             vertices.size(),
                             GL::Buffer::MapFlag::Write);

    CORRADE_CONSTEXPR_ASSERT(data, "could not map vertex data");
    Utility::copy(phasefieldData.meshData.vertexData(), data);
    vertices.unmap();
}

void UpdateScene::reuploadIndices()
{
    Debug{} << "Reuploading Indices";
    GL::Buffer& indices = phasefieldData.indices;
    auto data = indices.map(0,
                            indices.size(),
                            GL::Buffer::MapFlag::Write);

    CORRADE_CONSTEXPR_ASSERT(data, "could not map vertex data");
    Utility::copy(phasefieldData.meshData.indexData(), data);
    indices.unmap();
}

void UpdateScene::reupload(){
    reuploadVertices();
    reuploadIndices();
}