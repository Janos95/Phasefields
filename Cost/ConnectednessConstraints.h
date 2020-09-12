//
// Created by janos on 21.05.20.
//

#pragma once

#include "Functional.h"
#include "Dijkstra.h"
#include "StoppingCriteria.h"

namespace Phasefield {

namespace Mn = Magnum;
namespace Cr = Corrade;

struct ConnectednessConstraint {
    ConnectednessConstraint(Mesh&);

    template<class Scalar>
    void operator()(ArrayView<const Scalar> const& parameters,
                    ArrayView<const Scalar> const& weights,
                    Scalar& out,
                    ArrayView<Scalar> const& gradP,
                    ArrayView<Scalar> const& gradW);

    uint32_t numParameters() const;

    uint32_t numResiduals() const;

    struct Edge {
        Edge(Mg::UnsignedInt v1_, Mg::UnsignedInt v2_) : v1(std::min(v1_, v2_)), v2(std::max(v1_, v2_)) {}

        Mg::UnsignedInt v1, v2;

        auto operator<=>(const Edge&) const = default;
    };

    Array<const UnsignedInt> const& indices;
    Array<const Vector3d> const& vertices;
    Array<Edge> dualEdges;

    Array<double> lineElements;
    Array<double> areas;
    Array<double> diams;

    Array<Dijkstra<Face, Scalar> dijkstras;
    Array<StoppingCriteria> stops;
    Array<size_t> roots;
    size_t numComponents;
    Array<size_t> components;

    /* this is called holding our 'global' mutex from the gui thread */
    void updateVisualization(VisualizationProxy&);

    void drawImGuiOptions(bool&, DrawableType&, bool&);

    /* collects data for the updateVisualization callback
     * Always called after a call to evaluate.
     */
    //void collectVisualizationData(Containers::ArrayView<const double*> const& grad);

    double a = 0.05, b = 1.;
    double pathThickness = 0.01;

    //Cr::Containers::Array<InstanceData> instanceData; //tf from cylinder to path section
    //Paths* paths = nullptr;
    //bool updateInstanceData = false;
    //bool generateLineStrips = false;
    //bool updateGrad = false;
    //bool updateComponents = false;
    //bool updateWs = false;
};


extern template void ConnectednessConstraint::operator()(ArrayView<const double> const& parameters,
                                                 ArrayView<const double> const& weights,
                                                 double& out,
                                                 ArrayView<double> const& gradP,
                                                 ArrayView<double> const& gradW);

extern template Functional::Functional(ConnectednessConstraint);

}