//shamelessly stolen from https://github.com/mosra/magnum

#pragma once

#include <Magnum/GL/AbstractShaderProgram.h>
#include <Magnum/Shaders/Generic.h>
#include <Magnum/Shaders/visibility.h>

namespace Mg = Magnum;

namespace shaders {

class Phong : public Mg::GL::AbstractShaderProgram {
public:

    typedef Magnum::Shaders::Generic3D::Position Position;

    typedef Magnum::Shaders::Generic3D::Normal Normal;

    typedef Magnum::Shaders::Generic3D::Tangent Tangent;

    typedef Magnum::Shaders::Generic3D::TextureCoordinates TextureCoordinates;

    typedef Magnum::Shaders::Generic3D::Color3 Color3;

    typedef Magnum::Shaders::Generic3D::Color4 Color4;

    typedef Magnum::Shaders::Generic3D::TransformationMatrix TransformationMatrix;

    typedef Magnum::Shaders::Generic3D::NormalMatrix NormalMatrix;

    typedef typename Magnum::Shaders::Generic3D::TextureOffset TextureOffset;

    enum : Mg::UnsignedInt {
        ColorOutput = Magnum::Shaders::Generic3D::ColorOutput,
    };

    enum class Flag : Mg::UnsignedShort {

        AmbientTexture = 1 << 0,

        DiffuseTexture = 1 << 1,

        SpecularTexture = 1 << 2,

        NormalTexture = 1 << 4,

        AlphaMask = 1 << 3,

        VertexColor = 1 << 5,

        TextureTransformation = 1 << 6,

        InstancedTransformation = 1 << 9,

        InstancedTextureOffset = (1 << 10) | TextureTransformation
    };


    typedef Mg::Containers::EnumSet<Flag> Flags;

    explicit Phong(Flags flags = {}, Mg::UnsignedInt lightCount = 1);

    explicit Phong(Mg::NoCreateT) noexcept: Mg::GL::AbstractShaderProgram{Mg::NoCreate} {}

    /** @brief Copying is not allowed */
    Phong(const Phong&) = delete;

    /** @brief Move constructor */
    Phong(Phong&&) noexcept = default;

    /** @brief Copying is not allowed */
    Phong& operator=(const Phong&) = delete;

    /** @brief Move assignment */
    Phong& operator=(Phong&&) noexcept = default;

    /** @brief Flags */
    Flags flags() const { return _flags; }

    /** @brief Light count */
    Mg::UnsignedInt lightCount() const { return _lightCount; }

    Phong& setAmbientColor(const Magnum::Color4& color);

    Phong& bindAmbientTexture(Mg::GL::Texture2D& texture);

    Phong& setDiffuseColor(const Magnum::Color4& color);

    Phong& bindDiffuseTexture(Mg::GL::Texture2D& texture);

    Phong& bindNormalTexture(Mg::GL::Texture2D& texture);

    Phong& setSpecularColor(const Magnum::Color4& color);

    Phong& bindSpecularTexture(Mg::GL::Texture2D& texture);

    Phong& bindTextures(Mg::GL::Texture2D* ambient, Mg::GL::Texture2D* diffuse, Mg::GL::Texture2D* specular,
                        Mg::GL::Texture2D* normal
#ifdef MAGNUM_BUILD_DEPRECATED
                        = nullptr
#endif
    );

    Phong& setShininess(Mg::Float shininess);

    Phong& setAlphaMask(Mg::Float mask);

    Phong& setTransformationMatrix(const Mg::Matrix4& matrix);

    Phong& setNormalMatrix(const Mg::Matrix3x3& matrix);

    Phong& setProjectionMatrix(const Mg::Matrix4& matrix);

    Phong& setTextureMatrix(const Mg::Matrix3& matrix);

    Phong& setLightPositions(Mg::Containers::ArrayView<const Mg::Vector3> lights);

    Phong& setLightPositions(std::initializer_list<Mg::Vector3> lights);

    Phong& setLightPosition(Mg::UnsignedInt id, const Mg::Vector3& position);

    Phong& setLightPosition(const Mg::Vector3& position) {
        return setLightPositions({&position, 1});
    }

    Phong& setLightColors(Mg::Containers::ArrayView<const Magnum::Color4> colors);

    Phong& setLightColors(std::initializer_list<Magnum::Color4> colors);

    Phong& setLightColor(Mg::UnsignedInt id, const Magnum::Color4& color);

    Phong& setLightColor(const Magnum::Color4& color) {
        return setLightColors({&color, 1});
    }

private:

    Flags _flags;
    Mg::UnsignedInt _lightCount;
    Mg::Int _transformationMatrixUniform{0},
            _projectionMatrixUniform{1},
            _normalMatrixUniform{2},
            _textureMatrixUniform{3},
            _ambientColorUniform{4},
            _diffuseColorUniform{5},
            _specularColorUniform{6},
            _shininessUniform{7},
            _alphaMaskUniform{8};
    Mg::Int _lightPositionsUniform{10},
            _lightColorsUniform; /* 10 + lightCount, set in the constructor */
};

CORRADE_ENUMSET_OPERATORS(Phong::Flags)

}
