//
// Created by janos on 21.05.20.
//

#pragma once

#include "Functional.h"
#include "Dijkstra.h"
#include "Surface.h"

namespace Phasefield {

namespace Mn = Magnum;
namespace Cr = Corrade;

struct ConnectednessConstraint {
    explicit ConnectednessConstraint(Mesh&);

    template<class Scalar>
    void operator()(ArrayView<const Scalar> const& parameters,
                    ArrayView<const Scalar> const& weights,
                    Scalar& out,
                    ArrayView<Scalar> const& gradP,
                    ArrayView<Scalar> const& gradW);

    [[nodiscard]] size_t numParameters() const;

    Mesh* mesh;

    /* this is called holding our 'global' mutex from the gui thread */
    //void updateVisualization(VisualizationProxy&);

    //void drawImGuiOptions(bool&, DrawableType&, bool&);

    /* collects data for the updateVisualization callback
     * Always called after a call to evaluate.
     */
    //void collectVisualizationData(Containers::ArrayView<const double*> const& grad);

    double a = 0.05, b = 1.;
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

}