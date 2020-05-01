//
// Created by janos on 02.04.20.
//

#pragma once

#include "drawable_data.hpp"
#include "drawables.hpp"
#include "functional.hpp"
#include "problem.hpp"

#include <Corrade/Containers/Array.h>
#include <Corrade/Containers/Pointer.h>

#include <Magnum/Magnum.h>
#include <Magnum/Math/Vector3.h>
#include <Magnum/MeshTools/Duplicate.h>
#include <Magnum/Trade/MeshData.h>

namespace Cr = Corrade;
namespace Mg = Magnum;

enum ShaderType {
    FlatTextured = 0,
    PhongDiffuse = 1,
    MeshVisualizer = 2,
    MeshVisualizerObjectId = 3,
    MeshVisualizerTangent = 4,
    MeshVisualizerNormal = 5,
    MeshVisualizerFull = 6,
    VertexColor =7
};

struct PhasefieldData : Object3D, DrawableData {

    explicit PhasefieldData(Object3D * parent) : Object3D(parent) {}

    enum class Status: Magnum::UnsignedShort {
        NothingChanged,
        PrimitiveNewMesh,
        Subdivided,
        OptimizationUpdated,
        BrushUpdated
    };

    Cr::Containers::Array<Magnum::Vector3d> V;
    Cr::Containers::Array<Magnum::UnsignedInt> F;
    Cr::Containers::Array<Magnum::Double> phasefield;

    Mg::Trade::MeshData original{Magnum::MeshPrimitive::Points, 0};
    Mg::Trade::MeshData meshData{Magnum::MeshPrimitive::Points, 0};
    Mg::UnsignedInt numSubdivisions = 0;

    Status status = Status::NothingChanged;

    Cr::Containers::Array<Cr::Containers::Pointer<Mg::GL::AbstractShaderProgram>> shaders;
    Cr::Containers::Array<Magnum::GL::Texture2D> textures;

    solver::Problem problem;

    //void serialize(std::string const& path){
    //    auto sizeOriginal = original.serializedSize();
    //    auto sizeMeshData = meshData.serializedSize();

    //    Magnum::Containers::Array<char> buffer(Magnum::Containers::NoInit, sizeOriginal + sizeMeshData + sizeof(int));

    //    auto current = buffer.begin();
    //    original.serializeInto({current, sizeOriginal});
    //    current += sizeOriginal;
    //    meshData.serializeInto({current, sizeMeshData});
    //    std::memcpy(current, &numSubdivisions, sizeof(int));

    //    auto* file = std::fopen(path.c_str(), "wb");
    //    std::fwrite(buffer.data(), sizeof(char), buffer.size(), file);
    //    std::fclose(file);
    //}

    //void deserialize(std::string const& path){
    //    auto* file = std::fopen(path.c_str(), "rb");
    //    std::fseek(file, 0, SEEK_END); // seek to end
    //    auto filesize = std::ftell(file);

    //    std::fseek(file, 0, SEEK_SET); // seek to start
    //    Magnum::Containers::Array<char> buffer(Magnum::Containers::NoInit, filesize);
    //    std::fread(buffer.data(), sizeof(uint8_t), buffer.size(), file);

    //    Magnum::Trade::MeshData* p;

    //    status = Status::NewMesh;
    //}
};