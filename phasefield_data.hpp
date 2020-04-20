//
// Created by janos on 02.04.20.
//

#pragma once

#include "drawable_data.hpp"
#include "drawables.hpp"
#include "functional.hpp"

#include <Corrade/Containers/Array.h>
#include <Corrade/Containers/Pointer.h>

#include <Magnum/Magnum.h>
#include <Magnum/Math/Vector3.h>
#include <Magnum/MeshTools/Duplicate.h>

namespace Cr = Corrade;
namespace Mg = Magnum;

enum FunctionalType {
    DirichletEnergy = 0,
    DoubleWellPotential = 1,
    Area = 2,
    Connectedness = 3
};

enum ShaderType {
    FlatTextured = 0,
    PhongDiffuse = 1,
    MeshVisualizer = 2,
    VertexColor = 3
};

struct FunctionalParameters { char const* name; Mg::Float weight; };

struct PhasefieldData : ColorMapDrawableData {

    enum class Status: Magnum::UnsignedShort {
        NothingChanged,
        NewMesh,
        Subdivided,
        PhasefieldUpdated
    };

    Cr::Containers::Array<Magnum::Vector3> V;
    Cr::Containers::Array<Magnum::UnsignedInt> F;
    Cr::Containers::Array<Magnum::Float> phasefield;

    Mg::Trade::MeshData original{Magnum::MeshPrimitive::Points, 0};
    Mg::Trade::MeshData meshData{Magnum::MeshPrimitive::Points, 0};
    Mg::UnsignedInt numSubdivisions = 0;
    Status status = Status::NothingChanged;

    Drawable* drawable = nullptr;
    DrawableType type;

    struct Shader { Cr::Containers::Pointer<Mg::GL::AbstractShaderProgram> shader; ShaderType type;};
    Cr::Containers::Array<Shader> shaders;

    Cr::Containers::Array<Cr::Containers::Pointer<Functional>> functionals;

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