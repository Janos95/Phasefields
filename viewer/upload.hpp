//
// Created by janos on 13.03.20.
//

#include <Magnum/Trade/MeshData.h>
#include <Magnum/ImageView.h>

#include "drawable_data.hpp"

enum class CompileFlag: Magnum::UnsignedShort {
    GenerateFlatNormals = 1u << 0u,
    GenerateSmoothNormals = 1u << 1u,
    AddColorAttribute = 1u << 2u,
    AddNormalAttribute = 1u << 3u,
    AddTextureCoordinates = 1u << 4u
};

using CompileFlags = Corrade::Containers::EnumSet<CompileFlag>;

CORRADE_ENUMSET_OPERATORS(CompileFlags)

Magnum::Trade::MeshData preprocess(Magnum::Trade::MeshData const& meshData, CompileFlags flags = {});

void upload(DrawableData& drawableData, Magnum::Trade::MeshData& meshData);