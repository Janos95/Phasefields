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

struct LoadPrimitives : Viewer::AbstractEventHandler {
    PhasefieldData& phasefieldData;
    bool track = false;

    explicit LoadPrimitives(PhasefieldData& data): phasefieldData(data) {}

    struct AbstractPrimitiveOptions { virtual ~AbstractPrimitiveOptions() = default; };

    struct ComboElement {
        std::string name;
        PrimitiveType type;
        std::unique_ptr<AbstractPrimitiveOptions> options;
    };

    struct CapsuleOptions : AbstractPrimitiveOptions {
        std::uint32_t hemisphereRings = 30; //	Number of (face) rings for each hemisphere. Must be larger or equal to 1.
        std::uint32_t cylinderRings = 30; // 	Number of (face) rings for cylinder. Must be larger or equal to 1.
        std::uint32_t segments = 30; //	Number of (face) segments. Must be larger or equal to 3.
        float radius = 1.f;
        float length = 2.f;
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



