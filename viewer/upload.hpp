//
// Created by janos on 13.03.20.
//

#include <Magnum/Trade/MeshData.h>
#include <Magnum/ImageView.h>

#include "drawable_data.hpp"

enum class CompileFlag: Magnum::UnsignedShort {
    GenerateFlatNormals = 1,
    GenerateSmoothNormals = 2,
    AddColorAttribute = 3,
    AddNormalAttribute = 4,
    AddVertexCoordsAttribute = 5,
};

using CompileFlags = Corrade::Containers::EnumSet<CompileFlag>;

CORRADE_ENUMSET_OPERATORS(CompileFlags)

Magnum::Trade::MeshData preprocess(Magnum::Trade::MeshData& meshData, CompileFlag flags = {});

void upload(DrawableData& drawableData);