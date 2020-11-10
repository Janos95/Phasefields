//
// Created by janos on 21.05.20.
//

#pragma once

#include "Functional.h"
#include "Surface.h"

class adouble;

namespace Phasefield {

namespace Mn = Magnum;
namespace Cr = Corrade;

struct Node;

struct ConnectednessConstraint {
    explicit ConnectednessConstraint(Mesh&);

    template<class Scalar>
    void operator()(ArrayView<const Scalar> parameters,
                    ArrayView<const Scalar> weights,
                    Scalar& out,
                    ArrayView<Scalar> gradP,
                    ArrayView<Scalar> gradW);

    [[nodiscard]] size_t numParameters() const;
    [[nodiscard]] static FunctionalType::Value type() { return FunctionalType::ConnectednessConstraint; }

    void drawImGuiOptions(VisualizationProxy&);

    void draw(Node&);
    char const* getFormattedInterval();
    void recalculateInterval();

    void saveParameters(Cr::Utility::ConfigurationGroup&) const;
    void loadParameters(Cr::Utility::ConfigurationGroup const&);

    Mesh* mesh;

    int purePhase = 1;
    double edge0 = 0, edge1 = 1;
    double s = 0.5;
    bool drawComponents = false;
    bool drawGradient = false;
    bool ignoreSmallComponents = false;
    double* epsilon = nullptr;

    //double pathThickness = 0.01;

    //Cr::Containers::Array<InstanceData> instanceData; //tf from cylinder to path section
    //Paths* paths = nullptr;
    //bool updateInstanceData = false;
    //bool generateLineStrips = false;
    //bool updateGrad = false;
    //bool updateComponents = false;
    //bool updateWs = false;
};


DECLARE_FUNCTIONAL_CONSTRUCTOR(ConnectednessConstraint)
DECLARE_FUNCTIONAL_OPERATOR(ConnectednessConstraint, double)

#ifdef PHASEFIELD_WITH_ADOLC
DECLARE_FUNCTIONAL_OPERATOR(ConnectednessConstraint, adouble)
#endif

}

