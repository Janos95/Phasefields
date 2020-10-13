//
// Created by janos on 10/9/20.
//

#pragma once

#include "Types.h"

#include <Magnum/GL/AbstractShaderProgram.h>

namespace Phasefield {

namespace Mg = Magnum;

class DepthReinterpretShader : public GL::AbstractShaderProgram {
public:
    explicit DepthReinterpretShader(Mg::NoCreateT) : GL::AbstractShaderProgram{Mg::NoCreate} {}

    explicit DepthReinterpretShader();

    DepthReinterpretShader& bindDepthTexture(GL::Texture2D& texture);

};

}
