//
// Created by janos on 10/9/20.
//

#include "DepthReinterpretShader.h"

#include <Corrade/Utility/Resource.h>
#include <Corrade/Containers/Reference.h>
#include <Magnum/GL/Version.h>
#include <Magnum/GL/Shader.h>
#include <Magnum/GL/Texture.h>

namespace Cr = Corrade;
namespace Mg = Magnum;

namespace Phasefield {

DepthReinterpretShader::DepthReinterpretShader() {
    GL::Shader vert{GL::Version::GLES300, GL::Shader::Type::Vertex};
    GL::Shader frag{GL::Version::GLES300, GL::Shader::Type::Fragment};

    Cr::Utility::Resource rs{"data"};
    vert.addSource(rs.get("DepthReinterpretShader.vert"));
    frag.addSource(rs.get("DepthReinterpretShader.frag"));

    CORRADE_INTERNAL_ASSERT_OUTPUT(Mg::GL::Shader::compile({vert, frag}));

    attachShaders({vert, frag});
    CORRADE_INTERNAL_ASSERT_OUTPUT(link());

    setUniform(uniformLocation("depthTexture"), 7);
}

DepthReinterpretShader& DepthReinterpretShader::bindDepthTexture(GL::Texture2D& texture) {
    texture.bind(7);
    return *this;
}

}