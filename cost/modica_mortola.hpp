//
// Created by janos on 27.11.19.
//

#pragma once

#include "functional.h"
#include "c1_functions.hpp"
#include "stopping_criteria.hpp"
#include "fem.hpp"
#include "meta_data.hpp"
#include "permute.hpp"
#include "strang_rules.hpp"
#include "y_combinator.hpp"
#include "phasefield_tree.hpp"
#include "types.hpp"
#include "normalizeInto.hpp"

#include <Corrade/Containers/ArrayView.h>
#include <Corrade/Containers/GrowableArray.h>
#include <Corrade/Containers/StaticArray.h>
#include <Corrade/Utility/Algorithms.h>

#include <Magnum/Trade/MeshData.h>
#include <Magnum/Math/FunctionsBatch.h>

#include <Eigen/SparseCore>

#include <mutex>
#include <numeric>
#include <algorithm>

struct DirichletEnergy
{
    DirichletEnergy(
            Containers::Array<const Mg::Vector3d> const& vertices,
            Containers::Array<const Mg::UnsignedInt> const& indices);

    uint32_t numParameters() const ;

    template<class Scalar>
    bool evaluate(double const* parameters,
                  double* residual,
                  double** jacobians) const ;

    auto triangles() { return arrayCast<Mg::Vector3ui>(indices); }

    Containers::Array<const Mg::UnsignedInt> const& indices;
    Containers::Array<const Mg::Vector3d> const& vertices;
    Containers::Array<Mg::Double> areas;

    Eigen::SparseMatrix<double> stiffnessMatrix;
    SparseMatrix stiffnessMatrix2;
    Eigen::Matrix<double, Eigen::Dynamic, 1> diagonal;
};

template<template <class> class F>
struct IntegralFunctional
{
    IntegralFunctional(
            Containers::Array<const Mg::Vector3d> const& vs,
            Containers::Array<const Mg::UnsignedInt> const& is):
                Functional(GradientMetaData::fromLossAndName(L{}, "Integral Functional")),
                indices(is),
                vertices(vs),
                integralOperator(computeIntegralOperator(triangles(), vertices))
    {
    }

    uint32_t numParameters() const  {
        return vertices.size();
    }

    void updateIntegralOperator()  {
        integralOperator = computeIntegralOperator(triangles(), vertices);
    }

    template<class Scalar>
    bool evaluate(Scalar const* params,
                  Scalar* cost,
                  SparseMatrix* jacobians,
                  SparseMatrix* hessian) const {

        F<Scalar> f;
        *cost = 0.;
        for (int i = 0; i < vertices.size(); ++i){
            *cost += f.eval(params[i]) * integralOperator[i];
            if(jacobians)
                jacobians->values[i] = f.grad(params[i]) * integralOperator[i];
            if(hessian)
                hessian->values[i] = f.hess(params[i]) * integralOperator[i];
        }
        return true;
    }

    auto triangles() { return arrayCast<Mg::Vector3ui>(indices); }

    Containers::Array<const Mg::UnsignedInt> const& indices;
    Containers::Array<const Mg::Vector3d> const& vertices;
    Containers::Array<Mg::Double> integralOperator;
};


template<template <class> class F, class L>
struct AreaRegularizer final : IntegralFunctional<F> {

    AreaRegularizer(
            Containers::Array<Magnum::Vector3d> vs,
            Containers::Array<const Mg::UnsignedInt> const& is) :
        IntegralFunctional<F>(vs, is)
    {
        auto areas = computeAreas(triangles(), this->vertices);
        area = std::accumulate(areas.begin(), areas.end(), 0.);
    }

    AreaMetaData& getMeta(){ return *dynamic_cast<AreaMetaData*>(this->meta); }

    void updateInternalDataStructures()  {
        IntegralFunctional<F>::updateInternalDataStructures();
        auto areas = computeAreas(triangles(), this->vertices);
        area = std::accumulate(areas.begin(), areas.end(), 0.);
    }

    template<class Scalar>
    bool evaluate(Scalar const* params,
                  Scalar* cost,
                  double* jacobians) const {
        IntegralFunctional<Scalar, F, L>::evaluate(params, &currentArea, jacobians);
        *cost = currentArea - getMeta().areaRatio * area; /*@todo racy */
        return true;
    }

    auto triangles() { return arrayCast<Mg::Vector3ui>(this->indices); }

    mutable Mg::Double currentArea;
    Mg::Double area;
};


template<class Scalar, template<class> class SmoothStepFunc, class L>
struct AreaConstraints final : Functional {

    AreaMetaData& getMeta(){ return *dynamic_cast<AreaMetaData*>(this->meta); }

    AreaConstraints(Cr::Containers::ArrayView<const Mg::Vector3d> const& vs,
                    Cr::Containers::ArrayView<const Mg::Vector3ui> const& ts,
                    PhasefieldTree& t) : tree(t), integralOperator(computeIntegralOperator(ts,vs))
    {
    }

    Mg::UnsignedInt numParameters() const  { return tree.phasefieldData.size(); }
    Mg::UnsignedInt numResiduals() const  { return tree.numLeafs * 2; }


    bool evaluate(Scalar const* params,
                  Scalar* cost,
                  Scalar* jacobians) const {

        auto numParams = numParameters();
        auto numRes = numResiduals();
        auto size = tree.phasefieldSize;
        auto numPhasefields = tree.nodes.size();
        auto phasefieldSize = tree.phasefieldSize;
        auto& nodes = tree.nodes;

        DoubleView2D phasefields{{params, numParams}, {numPhasefields, phasefieldSize}};
        //DoubleView3D jac{jacobians, {numRes, numPhasefields, phasefieldSize}}; /* column major */

        auto visitor = YCombinator{
            [&](auto&& visitor, PhasefieldNode& node, Containers::Array<Scalar>& p) -> void {
                auto depth = node.depth;
                auto idx = node.idx;
                bool hasChildren = node.leftChild != PhasefieldNode::None;

                Containers::Array<Scalar> p1(hasChildren ? phasefieldSize : 0);
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

        Containers::Array<Scalar> p(phasefieldSize);
        std::fill(p.begin(), p.end(), 1.);
        std::fill_n(cost, numRes, 0.);
        visitor(tree.root(), p);

        return true;
    }

    SmoothStepFunc<Scalar> smoothStepFunc;

    PhasefieldTree& tree;
    Containers::Array<const Mg::Vector3ui> const& triangles;
    Containers::Array<const Mg::Vector3d> const& vertices;
    Containers::Array<Mg::Double> integralOperator;
};

using AreaRegularizer1 = AreaRegularizer<Mg::Double, Indicator, QuadraticLoss>;
using AreaRegularizer2 = AreaRegularizer<Mg::Double, SmootherStep, QuadraticLoss>;
using MultipleAreaConstraints = AreaConstraints<Mg::Double, SmootherStep, QuadraticLoss>;
using DoubleWellPotential = IntegralFunctional<Mg::Double, DoubleWell, TrivialLoss>;
