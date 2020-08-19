//
// Created by janos on 13.03.20.
//

#pragma once

#include <Magnum/Trade/Trade.h>
#include <Magnum/GL/GL.h>
#include <Corrade/Containers/EnumSet.h>
#include <Corrade/Containers/Containers.h>

namespace Mg = Magnum;
namespace Cr = Corrade;

enum class CompileFlag : Mg::UnsignedInt {
    GenerateFlatNormals = 1u << 0u,
    GenerateSmoothNormals = 1u << 1u,
    AddColorAttribute = 1u << 2u,
    AddNormalAttribute = 1u << 3u,
    AddTextureCoordinates = 1u << 4u,
};

using CompileFlags = Cr::Containers::EnumSet<CompileFlag>;

CORRADE_ENUMSET_OPERATORS(CompileFlags)

Mg::Trade::MeshData preprocess(
        Cr::Containers::ArrayView<const Mg::Vector3d> const& vertices,
        Cr::Containers::ArrayView<const Mg::UnsignedInt> const& indices,
        CompileFlags flags = {});

void upload(Mg::GL::Mesh& mesh, Mg::GL::Buffer& vertices, Mg::GL::Buffer& indices, Mg::Trade::MeshData& meshData);

void reuploadVertices(Mg::GL::Buffer& vertices, Mg::Trade::MeshData const& meshData);

