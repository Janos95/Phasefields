//
// Created by janos on 08.05.20.
//

#include "Paths.h"

#include <Corrade/Containers/GrowableArray.h>

#include <Magnum/MeshTools/Compile.h>
#include <Magnum/Primitives/Cylinder.h>
#include <Magnum/Trade/MeshData.h>

Paths::Paths():
        shader{Mg::Shaders::Phong::Flag::VertexColor|Mg::Shaders::Phong::Flag::InstancedTransformation},
        cylinder{Mg::MeshTools::compile(Mg::Primitives::cylinderSolid(3,5,0.5f))}
{
    shader.setAmbientColor(0x111111_rgbf)
            .setSpecularColor(0x330000_rgbf)
            .setLightPosition({10.0f, 15.0f, 5.0f});

    /* cylinder mesh, with an (initially empty) instance buffer */
    cylinder.addVertexBufferInstanced(instanceBuffer, 1, 0,
                                      Mg::Shaders::Phong::TransformationMatrix{},
                                      Mg::Shaders::Phong::NormalMatrix{},
                                      Mg::Shaders::Phong::Color3{});

}

void Paths::draw(const Mg::Matrix4& transformation, Mg::Matrix4 const& projection){
    if(!drawPaths || instanceData.empty()) return;
    Mg::Containers::arrayResize(instanceDataTransformed, Mg::Containers::NoInit, instanceData.size());

    for (int i = 0; i < instanceData.size(); ++i) {
        instanceDataTransformed[i].normalMatrix = transformation.normalMatrix() * instanceData[i].normalMatrix;
        instanceDataTransformed[i].tf = transformation * instanceData[i].tf;
        instanceDataTransformed[i].color = instanceData[i].color;
    }

    instanceBuffer.setData(instanceDataTransformed, Mg::GL::BufferUsage::DynamicDraw);
    cylinder.setInstanceCount(instanceDataTransformed.size());
    shader.setProjectionMatrix(projection)
            .draw(cylinder);
}


