//
// Created by janos on 03.04.20.
//

#include "update_scene.hpp"

void UpdateScene::tickEvent(Scene& scene)  {

    if(!phasefieldData) return;

    switch (phasefieldData->status) {
        case PhasefieldData::Status::PhasefieldUpdated: reuploadVertices(*scene); break;
        case PhasefieldData::Status::NewMesh: uploadMesh(*scene); break;
        case PhasefieldData::Status::Subdivided: break;
    }

    scene->setDirty();

    bool newMeshData = true;
    if(m_newMeshData.compare_exchange_strong(newMeshData, false)) {
        Debug{} << "Uploading new mesh";
        if(scene) scene->reset();
        else scene = new Scene(); //@todo scene leaks, but I am not sure who should own it...
        {
            std::lock_guard l(m_mutex);
            scene->addObject("mesh", upload(*m_meshData));
        }
        scene->addNode("mesh", "mesh", ShaderType::Phong);
        scene->setDirty();
    }
    bool reupload = true;
    if(m_reupload.compare_exchange_strong(reupload, false)){
        auto* obj = scene->getObject("mesh");
        CORRADE_INTERNAL_ASSERT(obj);
        GL::Buffer& vertices = obj->vertices;
        auto data = vertices.map(0,
                                 vertices.size(),
                                 GL::Buffer::MapFlag::Write);

        CORRADE_CONSTEXPR_ASSERT(data, "could not map vertex data");
        {
            std::lock_guard l(m_mutex);
            Utility::copy(m_meshData->vertexData(), data);
        }
        vertices.unmap();
        scene->setDirty();
    }
}

void UpdateScene::makeScene(Scene& scene){
    scene.reset();
    scene.addObject("mesh", upload())
    scene.addNode("mesh", "mesh", ShaderType::Phong);
    scene.setDirty();
}

void UpdateScene::reuploadVertices(Scene& scene)
{
    auto* obj = scene.getObject("mesh");
    CORRADE_INTERNAL_ASSERT(obj);
    GL::Buffer& vertices = obj->vertices;
    auto data = vertices.map(0,
                             vertices.size(),
                             GL::Buffer::MapFlag::Write);

    CORRADE_CONSTEXPR_ASSERT(data, "could not map vertex data");
    Utility::copy(m_meshData->vertexData(), data);
    vertices.unmap();
    scene.setDirty();
}

void UpdateScene::uploadIndices(Scene& scene)
{
    Debug{} << "Uploading new mesh";
    if(scene) scene->reset();
    else scene = new Scene(); //@todo scene leaks, but I am not sure who should own it...
    {
        std::lock_guard l(m_mutex);
        scene->addObject("mesh", upload(*m_meshData));
    }
    scene->addNode("mesh", "mesh", ShaderType::Phong);
}