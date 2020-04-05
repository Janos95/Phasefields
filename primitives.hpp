//
// Created by janos on 04.04.20.
//

#pragma once

#include "phasefield_data.hpp"
#include "viewer.hpp"


enum class PrimitiveType: Magnum::UnsignedShort {
    Capsule,
    U
};

Magnum::Trade::MeshData uShapedSquare();

class LoadPrimitives : Viewer::AbstractEventHandler {
    PhasefieldData* m_phasefieldData = nullptr;

public:

    struct AbstractPrimitiveOptions { virtual ~AbstractPrimitiveOptions() = default; };

    struct ComboElement {
        std::string name;
        PrimitiveType type;
        std::unique_ptr<AbstractPrimitiveOptions> options;
    };

    struct CapsuleOptions : AbstractPrimitiveOptions {
        int hemisphereRings = 1; //	Number of (face) rings for each hemisphere. Must be larger or equal to 1.
        int cylinderRings = 1; // 	Number of (face) rings for cylinder. Must be larger or equal to 1.
        int segments = 3; //	Number of (face) segments. Must be larger or equal to 3.
        float halfLength = 1.; //	Half the length of cylinder part
    };

    struct UOptions : AbstractPrimitiveOptions {
        float height = 1.f;
        float width = 1.f;
        float innerWidth = .5f;
        float innerHeight = .5f;
    };

    void drawImGui() override;
    void load(ComboElement&);
};



