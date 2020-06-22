//
// Created by janos on 27.11.19.
//

#pragma once

#include "functional.hpp"
#include "c1_functions.hpp"
#include "stopping_criteria.hpp"
#include "fem.hpp"
#include "meta_data.hpp"
#include "permute.hpp"
#include "strang_rules.hpp"
#include "y_combinator.hpp"
#include "phasefield_tree.hpp"
#include "types.hpp"

#include <Corrade/Containers/ArrayView.h>
#include <Corrade/Containers/GrowableArray.h>
#include <Corrade/Containers/StaticArray.h>
#include <Corrade/Utility/Algorithms.h>

#include <Magnum/Trade/MeshData.h>
#include <Magnum/Math/FunctionsBatch.h>

#include <Eigen/SparseCore>

#include <normalizeInto.hpp>

#include <mutex>
#include <numeric>
#include <algorithm>

namespace Cr = Corrade;
namespace Mg = Magnum;

using DoubleView3D = Cr::Containers::StridedArrayView3D<double>;
using DoubleView2D = Cr::Containers::StridedArrayView2D<double>;
using DoubleView1D = Cr::Containers::StridedArrayView1D<double>;
using DoubleView = Cr::Containers::ArrayView<double>;

template<class Scalar>
struct DirichletEnergy final : Functional<Scalar>
{
    DirichletEnergy(
            Array<const Mg::Vector3d> const& vertices,
            Array<const Mg::UnsignedInt> const& indices);

    uint32_t numParameters() const override;

    bool evaluate(double const* parameters,
                  double* residual,
                  double** jacobians) const override;

    auto triangles() { return arrayCast<Mg::Vector3ui>(indices); }

    Array<const Mg::UnsignedInt> const& indices;
    Array<const Mg::Vector3d> const& vertices;
    Array<Mg::Double> areas;

    Eigen::SparseMatrix<Scalar> stiffnessMatrix;
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> diagonal;
};

template<class Scalar, template <class> class F, class L>
struct IntegralFunctional : Functional<Scalar>
{
    IntegralFunctional(
            Array<const Mg::Vector3d> const& vs,
            Array<const Mg::UnsignedInt> const& is,
            F<Scalar> f_ = {}):
                Functional<Scalar>(GradientMetaData::fromLossAndName(L{}, "Integral Functional")),
                indices(is),
                vertices(vs),
                f((F<Scalar>&&)f_),
                integralOperator(computeIntegralOperator(triangles(), vertices))
    {
    }

    uint32_t numParameters() const override {
        return vertices.size();
    }

    void updateIntegralOperator() override {
        integralOperator = computeIntegralOperator(triangles(), vertices);
    }

    virtual bool evaluate(double const* params,
                  double* cost,
                  double* jacobians) const override{
        *cost = 0.;
        for (int i = 0; i < vertices.size(); ++i){
            *cost += f.eval(params[i]) * integralOperator[i];
            if(jacobians) jacobians[i] = f.grad(params[i]) * integralOperator[i];
        }

        return true;
    }

    auto triangles() { return arrayCast<Mg::Vector3ui>(indices); }

    F<Scalar> f;
    Array<const Mg::UnsignedInt> const& indices;
    Array<const Mg::Vector3d> const& vertices;
    Array<Mg::Double> integralOperator;
};


template<class Scalar, template <class> class F, class L>
struct AreaRegularizer final : IntegralFunctional<Scalar, F, L> {

    AreaRegularizer(
            Array<Magnum::Vector3d> vs,
            Array<const Mg::UnsignedInt> const& is) :
        IntegralFunctional<Scalar, F, L>(vs, is)
    {
        auto areas = computeAreas(triangles(), this->vertices);
        area = std::accumulate(areas.begin(), areas.end(), 0.);
    }

    AreaMetaData& getMeta(){ return *dynamic_cast<AreaMetaData*>(this->meta); }

    void updateInternalDataStructures() override {
        IntegralFunctional<Scalar, F, L>::updateInternalDataStructures();
        auto areas = computeAreas(triangles(), this->vertices);
        area = std::accumulate(areas.begin(), areas.end(), 0.);
    }

    bool evaluate(double const* params,
                  double* cost,
                  double* jacobians) const override{
        IntegralFunctional<Scalar, F, L>::evaluate(params, &currentArea, jacobians);
        *cost = currentArea - getMeta().areaRatio * area; /*@todo racy */
        return true;
    }

    auto triangles() { return arrayCast<Mg::Vector3ui>(this->indices); }

    mutable Mg::Double currentArea;
    Mg::Double area;
};


template<class Scalar, template<class> class SmoothStepFunc, class L>
struct AreaConstraints final : Functional<Scalar> {

    AreaMetaData& getMeta(){ return *dynamic_cast<AreaMetaData*>(this->meta); }

    AreaConstraints(Cr::Containers::ArrayView<const Mg::Vector3d> const& vs,
                    Cr::Containers::ArrayView<const Mg::Vector3ui> const& ts,
                    PhasefieldTree& t) : tree(t), integralOperator(computeIntegralOperator(ts,vs))
    {
    }

    Mg::UnsignedInt numParameters() const override { return tree.phasefieldData.size(); }
    Mg::UnsignedInt numResiduals() const override { return tree.numLeafs * 2; }


    bool evaluate(Scalar const* params,
                  Scalar* cost,
                  Scalar* jacobians) const override{

        auto numParams = numParameters();
        auto numRes = numResiduals();
        auto size = tree.phasefieldSize;
        auto numPhasefields = tree.nodes.size();
        auto phasefieldSize = tree.phasefieldSize;
        auto& nodes = tree.nodes;

        DoubleView2D phasefields{{params, numParams}, {numPhasefields, phasefieldSize}};
        //DoubleView3D jac{jacobians, {numRes, numPhasefields, phasefieldSize}}; /* column major */

        auto visitor = YCombinator{
            [&](auto&& visitor, PhasefieldNode& node, Array<Scalar>& p) -> void {
                auto depth = node.depth;
                auto idx = node.idx;
                bool hasChildren = node.leftChild != PhasefieldNode::None;

                Array<Scalar> p1(hasChildren ? phasefieldSize : 0);
                if(hasChildren)
                    Cr::Utility::copy(p, p1);


                for (uint32_t i = 0; i < size; ++i) {
                    Scalar pos = smoothStepFunc.eval(phasefields[node.leftChild][i]);
                    Scalar neg = smoothStepFunc.eval(-phasefields[node.leftChild][i]);

                    if(hasChildren){
                        p[i] *= pos;
                        p1[i] *= neg;
                    } else {
                        cost[2*idx] += pos * p[i];
                        cost[2*idx + 1] += neg * p[i];
                    }
                }

                if(hasChildren){
                    visitor(nodes[node.leftChild], p);
                    visitor(nodes[node.rightChild], p1);
                }
            }
        };

        Array<Scalar> p(phasefieldSize);
        std::fill(p.begin(), p.end(), 1.);
        std::fill_n(cost, numRes, 0.);
        visitor(tree.root(), p);

        return true;
    }

    SmoothStepFunc<Scalar> smoothStepFunc;

    PhasefieldTree& tree;
    Array<const Mg::Vector3ui> const& triangles;
    Array<const Mg::Vector3d> const& vertices;
    Array<Mg::Double> integralOperator;
};

using AreaRegularizer1 = AreaRegularizer<Mg::Double, Indicator, QuadraticLoss>;
using AreaRegularizer2 = AreaRegularizer<Mg::Double, SmootherStep, QuadraticLoss>;
using MultipleAreaConstraints = AreaConstraints<Mg::Double, SmootherStep, QuadraticLoss>;
using DoubleWellPotential = IntegralFunctional<Mg::Double, DoubleWell, TrivialLoss>;
