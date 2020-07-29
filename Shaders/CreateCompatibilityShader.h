//shamelessly stolen from https://github.com/mosra/magnum

#pragma once

#include <Corrade/Utility/Resource.h>

#include "Magnum/GL/Context.h"
#include "Magnum/GL/Extensions.h"
#include "Magnum/GL/Shader.h"

namespace Mg = Magnum;
namespace Cr = Corrade;

namespace Implementation {

inline Mg::GL::Shader createCompatibilityShader(const Cr::Utility::Resource& rs, Mg::GL::Version version, Mg::GL::Shader::Type type) {
    Mg::GL::Shader shader(version, type);

    #ifndef MAGNUM_TARGET_GLES
    if(Mg::GL::Context::current().isExtensionDisabled<Mg::GL::Extensions::ARB::explicit_attrib_location>(version))
        shader.addSource("#define DISABLE_GL_ARB_explicit_attrib_location\n");
    if(Mg::GL::Context::current().isExtensionDisabled<Mg::GL::Extensions::ARB::shading_language_420pack>(version))
        shader.addSource("#define DISABLE_GL_ARB_shading_language_420pack\n");
    if(Mg::GL::Context::current().isExtensionDisabled<Mg::GL::Extensions::ARB::explicit_uniform_location>(version))
        shader.addSource("#define DISABLE_GL_ARB_explicit_uniform_location\n");
    #endif

    #ifndef MAGNUM_TARGET_GLES2
    if(type == Mg::GL::Shader::Type::Vertex && Mg::GL::Context::current().isExtensionDisabled<Mg::GL::Extensions::MAGNUM::shader_vertex_id>(version))
        shader.addSource("#define DISABLE_GL_MAGNUM_shader_vertex_id\n");
    #endif

    /* My Android emulator (running on NVidia) doesn't define GL_ES
       preprocessor macro, thus *all* the stock shaders fail to compile */
    /** @todo remove this when Android emulator is sane */
    #ifdef CORRADE_TARGET_ANDROID
    shader.addSource("#ifndef GL_ES\n#define GL_ES 1\n#endif\n");
    #endif

    shader.addSource(rs.get("compatibility.glsl"));
    return shader;
}

}

