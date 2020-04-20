//
// Created by janos on 02.02.20.
//

#pragma once

#include <Corrade/Containers/GrowableArray.h>
#include <Corrade/Containers/Optional.h>
#include <Corrade/Containers/Pointer.h>

#include <Magnum/GL/Texture.h>
#include <Magnum/GL/Mesh.h>
#include <Magnum/Magnum.h>
#include <Magnum/Math/Color.h>
#include <Magnum/Trade/MeshData.h>
#include <Magnum/DebugTools/ColorMap.h>
#include <Magnum/PixelFormat.h>
#include <Magnum/ImageView.h>
#include <Magnum/GL/TextureFormat.h>
#include <Magnum/SceneGraph/Object.h>
#include <Magnum/SceneGraph/RigidMatrixTransformation3D.h>


using Object = Magnum::SceneGraph::Object<Magnum::SceneGraph::RigidMatrixTransformation3D>;

struct DrawableData : Object
{
    Magnum::GL::Buffer vertices, indices;
    Magnum::GL::Mesh mesh;

    virtual ~DrawableData() = default;
};

struct UniformColorDrawableData : virtual DrawableData
{
    Magnum::Color4 color;
};

enum class ColorMapType: Magnum::UnsignedShort {
    Turbo = 0,
    Magma = 1,
    Plasma = 2,
    Inferno = 3,
    Viridis = 4
};

struct ColorMapDrawableData : virtual DrawableData
{

    ColorMapDrawableData() {
        for(auto&& map : {
                            Magnum::DebugTools::ColorMap::turbo(),
                            Magnum::DebugTools::ColorMap::magma(),
                            Magnum::DebugTools::ColorMap::plasma(),
                            Magnum::DebugTools::ColorMap::inferno(),
                            Magnum::DebugTools::ColorMap::viridis()
                       })
        {
            const Magnum::Vector2i size{Magnum::Int(map.size()), 1};
            Magnum::GL::Texture2D texture;
            texture.setMinificationFilter(Magnum::SamplerFilter::Linear)
                    .setMagnificationFilter(Magnum::SamplerFilter::Linear)
                    .setWrapping(Magnum::SamplerWrapping::ClampToEdge) // or Repeat
                    .setStorage(1, Magnum::GL::TextureFormat::RGB8, size) // or SRGB8
                    .setSubImage(0, {}, Magnum::ImageView2D{Magnum::PixelFormat::RGB8Srgb, size, map});
            Corrade::Containers::arrayAppend(textures, std::move(texture));
        }
    }

    Corrade::Containers::Array<Magnum::GL::Texture2D> textures;
};






