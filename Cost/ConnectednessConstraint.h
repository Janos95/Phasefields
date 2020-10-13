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

    Mesh* mesh;

    /* this is called holding our 'global' mutex from the gui thread */
    //void updateVisualization(VisualizationProxy&);

    //void drawImGuiOptions(bool&, DrawableType&, bool&);

    /* collects data for the updateVisualization callback
     * Always called after a call to evaluate.
     */
    //void collectVisualizationData(Containers::ArrayView<const double*> const& grad);

    double a = -1.1, b = 0;
    bool drawComponents = false;
    bool drawGradient = false;
    bool ignoreSmallComponents = false;
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

