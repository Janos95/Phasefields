//shamelessly stolen from https://github.com/mosra/magnum

#include "Phong.h"

#ifdef MAGNUM_TARGET_GLES
#include <Corrade/Containers/Array.h>
#endif
#include <Corrade/Containers/EnumSet.hpp>
#include <Corrade/Containers/Reference.h>
#include <Corrade/Utility/FormatStl.h>
#include <Corrade/Utility/Resource.h>

#include <Magnum/GL/Context.h>
#include <Magnum/GL/Extensions.h>
#include <Magnum/GL/Shader.h>
#include <Magnum/GL/Texture.h>
#include <Magnum/Math/Color.h>
#include <Magnum/Math/Matrix3.h>
#include <Magnum/Math/Matrix4.h>

#include "CreateCompatibilityShader.h"

using namespace Magnum;
using namespace Corrade;

namespace Phasefield::Shaders {

namespace {
enum : Int {
    AmbientTextureUnit = 0,
    DiffuseTextureUnit = 1,
    SpecularTextureUnit = 2,
    NormalTextureUnit = 3
};
}

Phong::Phong(const Flags flags, const UnsignedInt lightCount) : _flags{flags}, _lightCount{lightCount},
                                                                _lightColorsUniform{
                                                                        _lightPositionsUniform + Int(lightCount)} {
    CORRADE_ASSERT(!(flags & Flag::TextureTransformation) || (flags & (Flag::AmbientTexture | Flag::DiffuseTexture |
                                                                       Flag::SpecularTexture | Flag::NormalTexture)),
                   "Shaders::Phong: texture transformation enabled but the shader is not textured",);

#ifdef MAGNUM_BUILD_STATIC
    /* Import resources on static build, if not already */
    if(!Utility::Resource::hasGroup("MagnumShaders"))
        importShaderResources();
#endif
    Utility::Resource rs("MagnumShaders");

#ifndef MAGNUM_TARGET_GLES
    const GL::Version version = GL::Context::current().supportedVersion(
            {GL::Version::GL320, GL::Version::GL310, GL::Version::GL300, GL::Version::GL210});
#else
    const GL::Version version = GL::Context::current().supportedVersion({GL::Version::GLES300, GL::Version::GLES200});
#endif

    GL::Shader vert = Implementation::createCompatibilityShader(rs, version, GL::Shader::Type::Vertex);
    GL::Shader frag = Implementation::createCompatibilityShader(rs, version, GL::Shader::Type::Fragment);

#ifndef MAGNUM_TARGET_GLES
    std::string lightInitializer;
    if(lightCount) {
        /* Initializer for the light color array -- we need a list of vec4(1.0)
           joined by commas. For GLES we'll simply upload the values directly. */
        constexpr const char lightInitializerPreamble[] = "#define LIGHT_COLOR_INITIALIZER ";
        constexpr std::size_t lightInitializerPreambleSize =
                Containers::arraySize(lightInitializerPreamble) - 1;
        constexpr const char lightInitializerItem[] = "vec4(1.0), ";
        constexpr std::size_t lightInitializerItemSize =
                Containers::arraySize(lightInitializerItem) - 1;
        lightInitializer.reserve(
                Containers::arraySize(lightInitializerPreamble) - 1 + lightCount*lightInitializerItemSize);
        lightInitializer.append(lightInitializerPreamble, lightInitializerPreambleSize);
        for(std::size_t i = 0; i != lightCount; ++i)
            lightInitializer.append(lightInitializerItem, lightInitializerItemSize);

        /* Drop the last comma and add a newline at the end */
        lightInitializer[lightInitializer.size() - 2] = '\n';
        lightInitializer.resize(lightInitializer.size() - 1);
    }
#endif

    vert.addSource(flags & (Flag::AmbientTexture | Flag::DiffuseTexture | Flag::SpecularTexture | Flag::NormalTexture)
                   ? "#define TEXTURED\n" : "")
        .addSource(flags & Flag::NormalTexture ? "#define NORMAL_TEXTURE\n" : "")
        .addSource(flags & Flag::VertexColor ? "#define VERTEX_COLOR\n" : "")
        .addSource(flags & Flag::TextureTransformation ? "#define TEXTURE_TRANSFORMATION\n" : "")
        .addSource(Utility::formatString("#define LIGHT_COUNT {}\n", lightCount))
        .addSource(flags & Flag::InstancedTransformation ? "#define INSTANCED_TRANSFORMATION\n" : "")
        .addSource(flags >= Flag::InstancedTextureOffset ? "#define INSTANCED_TEXTURE_OFFSET\n" : "")
        .addSource(rs.get("generic.glsl"))
        .addSource(rs.get("Phong.vert"));
    frag.addSource(flags & Flag::AmbientTexture ? "#define AMBIENT_TEXTURE\n" : "")
        .addSource(flags & Flag::DiffuseTexture ? "#define DIFFUSE_TEXTURE\n" : "")
        .addSource(flags & Flag::SpecularTexture ? "#define SPECULAR_TEXTURE\n" : "")
        .addSource(flags & Flag::NormalTexture ? "#define NORMAL_TEXTURE\n" : "")
        .addSource(flags & Flag::VertexColor ? "#define VERTEX_COLOR\n" : "")
        .addSource(Utility::formatString(
                "#define LIGHT_COUNT {}\n"
                "#define LIGHT_COLORS_LOCATION {}\n", lightCount, _lightPositionsUniform + lightCount));
#ifndef MAGNUM_TARGET_GLES
    if(lightCount) frag.addSource(std::move(lightInitializer));
#endif
    frag.addSource(rs.get("generic.glsl"))
        .addSource(rs.get("Phong.frag"));

    CORRADE_INTERNAL_ASSERT_OUTPUT(GL::Shader::compile({vert, frag}));

    attachShaders({vert, frag});

    /* ES3 has this done in the shader directly and doesn't even provide
       bindFragmentDataLocation() */
#if !defined(MAGNUM_TARGET_GLES) || defined(MAGNUM_TARGET_GLES2)
#ifndef MAGNUM_TARGET_GLES
    if(!GL::Context::current().isExtensionSupported<GL::Extensions::ARB::explicit_attrib_location>(version))
#endif
    {
        bindAttributeLocation(Position::Location, "position");
        if(lightCount)
            bindAttributeLocation(Normal::Location, "normal");
        if((flags & Flag::NormalTexture) && lightCount)
            bindAttributeLocation(Tangent::Location, "tangent");
        if(flags & Flag::VertexColor)
            bindAttributeLocation(Color3::Location, "vertexColor"); /* Color4 is the same */
        if(flags & (Flag::AmbientTexture | Flag::DiffuseTexture | Flag::SpecularTexture))
            bindAttributeLocation(TextureCoordinates::Location, "textureCoordinates");
        if(flags & Flag::InstancedTransformation)
            bindAttributeLocation(TransformationMatrix::Location, "instancedTransformationMatrix");
        if(flags >= Flag::InstancedTextureOffset)
            bindAttributeLocation(TextureOffset::Location, "instancedTextureOffset");
    }
#endif

    CORRADE_INTERNAL_ASSERT_OUTPUT(link());

#ifndef MAGNUM_TARGET_GLES
    if(!GL::Context::current().isExtensionSupported<GL::Extensions::ARB::explicit_uniform_location>(version))
#endif
    {
        _transformationMatrixUniform = uniformLocation("transformationMatrix");
        if(flags & Flag::TextureTransformation)
            _textureMatrixUniform = uniformLocation("textureMatrix");
        _projectionMatrixUniform = uniformLocation("projectionMatrix");
        _ambientColorUniform = uniformLocation("ambientColor");
        if(lightCount) {
            _normalMatrixUniform = uniformLocation("normalMatrix");
            _diffuseColorUniform = uniformLocation("diffuseColor");
            _specularColorUniform = uniformLocation("specularColor");
            _shininessUniform = uniformLocation("shininess");
            _lightPositionsUniform = uniformLocation("lightPositions");
            _lightColorsUniform = uniformLocation("lightColors");
        }
    }

#ifndef MAGNUM_TARGET_GLES
    if(flags && !GL::Context::current().isExtensionSupported<GL::Extensions::ARB::shading_language_420pack>(version))
#endif
    {
        if(flags & Flag::AmbientTexture) setUniform(uniformLocation("ambientTexture"), AmbientTextureUnit);
        if(lightCount) {
            if(flags & Flag::DiffuseTexture) setUniform(uniformLocation("diffuseTexture"), DiffuseTextureUnit);
            if(flags & Flag::SpecularTexture) setUniform(uniformLocation("specularTexture"), SpecularTextureUnit);
            if(flags & Flag::NormalTexture) setUniform(uniformLocation("normalTexture"), NormalTextureUnit);
        }
    }

    /* Set defaults in OpenGL ES (for desktop they are set in shader code itself) */
#ifdef MAGNUM_TARGET_GLES
    /* Default to fully opaque white so we can see the textures */
    if(flags & Flag::AmbientTexture) setAmbientColor(Magnum::Color4{1.0f});
    else setAmbientColor(Magnum::Color4{0.0f});
    setTransformationMatrix({});
    setProjectionMatrix({});
    if(lightCount) {
        setDiffuseColor(Magnum::Color4{1.0f});
        setSpecularColor(Magnum::Color4{1.0f});
        setShininess(80.0f);
        setLightColors(Containers::Array<Magnum::Color4>{Containers::DirectInit, lightCount, Magnum::Color4{1.0f}});
        /* Light position is zero by default */
        setNormalMatrix({});
    }
    if(flags & Flag::TextureTransformation) setTextureMatrix({});
    if(flags & Flag::AlphaMask) setAlphaMask(0.5f);
    /* Object ID is zero by default */
#endif
}



Phong& Phong::setShininess(Float shininess) {
    if(_lightCount) setUniform(_shininessUniform, shininess);
    return *this;
}

Phong& Phong::setTransformationMatrix(const Matrix4& matrix) {
    setUniform(_transformationMatrixUniform, matrix);
    return *this;
}

Phong& Phong::setNormalMatrix(const Matrix3x3& matrix) {
    if(_lightCount) setUniform(_normalMatrixUniform, matrix);
    return *this;
}

Phong& Phong::setProjectionMatrix(const Matrix4& matrix) {
    setUniform(_projectionMatrixUniform, matrix);
    return *this;
}

Phong& Phong::setLightPositions(const Containers::ArrayView<const Vector3> positions) {
    CORRADE_ASSERT(_lightCount == positions.size(),
                   "Shaders::Phong::setLightPositions(): expected" << _lightCount << "items but got"
                                                                   << positions.size(), *this);
    if(_lightCount) setUniform(_lightPositionsUniform, positions);
    return *this;
}

Phong& Phong::setLightPosition(UnsignedInt id, const Vector3& position) {
    CORRADE_ASSERT(id < _lightCount,
                   "Shaders::Phong::setLightPosition(): light ID" << id << "is out of bounds for" << _lightCount
                                                                  << "lights", *this);
    setUniform(_lightPositionsUniform + id, position);
    return *this;
}

/* It's light, but can't be in the header because MSVC needs to know the size
   of Vector3 for the initializer list use */
Phong& Phong::setLightPositions(std::initializer_list<Vector3> lights) {
    return setLightPositions({lights.begin(), lights.size()});
}

Phong& Phong::setLightColors(const Containers::ArrayView<const Magnum::Color4> colors) {
    CORRADE_ASSERT(_lightCount == colors.size(),
                   "Shaders::Phong::setLightColors(): expected" << _lightCount << "items but got" << colors.size(),
                   *this);
    if(_lightCount) setUniform(_lightColorsUniform, colors);
    return *this;
}

/* It's light, but can't be in the header because MSVC needs to know the size
   of Color for the initializer list use */
Phong& Phong::setLightColors(std::initializer_list<Magnum::Color4> colors) {
    return setLightColors({colors.begin(), colors.size()});
}

Phong& Phong::setLightColor(UnsignedInt id, const Magnum::Color4& color) {
    CORRADE_ASSERT(id < _lightCount,
                   "Shaders::Phong::setLightColor(): light ID" << id << "is out of bounds for" << _lightCount
                                                               << "lights", *this);
    setUniform(_lightColorsUniform + id, color);
    return *this;
}

}
