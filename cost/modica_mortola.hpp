//
// Created by janos on 27.11.19.
//

#pragma once

#include "functional.hpp"
#include "c1_functions.hpp"
#include "stopping_criteria.hpp"
#include "fem.hpp"

#include <Corrade/Containers/ArrayView.h>
#include <Corrade/Containers/GrowableArray.h>
#include <Magnum/Trade/MeshData.h>

#include <Eigen/SparseCore>
#include <numeric>

#include <normalizeInto.hpp>
#include <mutex>

namespace Cr = Corrade;
namespace Mg = Magnum;

struct GradientMetaData : Functional::MetaData {

    GradientMetaData(Mg::Trade::MeshData& md, VisualizationFlags& u, std::mutex& m) :
        meshData(md), update(u), mutex(m)
    {
    }

    void visualizeGradient(Cr::Containers::ArrayView<const Mg::Double> const& gradient) override {
        std::lock_guard l(mutex);
        if(flags & VisualizationFlag::Gradient){
            auto coords = meshData.mutableAttribute(Mg::Trade::MeshAttribute::TextureCoordinates);
            auto xcoords = Cr::Containers::arrayCast<2, Mg::Float>(coords).template slice<1>();
            normalizeInto(gradient, xcoords);
            update = VisualizationFlag::Gradient;
        }
    }

    Mg::Trade::MeshData& meshData;
    VisualizationFlags& update;
    std::mutex& mutex;
};


struct DirichletEnergy :  Functional
{
    DirichletEnergy(
            Cr::Containers::ArrayView<const Mg::Vector3d> const& vertices,
            Cr::Containers::ArrayView<const Mg::Vector3ui> const& faces);

    int numParameters() const override;

    bool evaluate(double const* parameters,
                  double* residual,
                  double* jacobians) const override;

    Cr::Containers::ArrayView<const Mg::Vector3ui> triangles;
    Cr::Containers::ArrayView<const Mg::Vector3d> vertices;
    Cr::Containers::Array<Mg::Double> areas;

    Eigen::SparseMatrix<Mg::Double> stiffnessMatrix;
    Eigen::VectorXd diagonal;
};

template<class F, class L, FunctionalType type>
struct IntegralFunctional : Functional
{
    IntegralFunctional(
            Cr::Containers::ArrayView<const Mg::Vector3d> const& vs,
            Cr::Containers::ArrayView<const Mg::Vector3ui> const& ts,
            F f_ = {}):
                Functional(Functional::MetaData::AllocateFromLoss(L{}), type),
                triangles(ts),
                vertices(vs),
                f((F&&)f_),
                integralOperator(computeIntegralOperator(ts,vs))
    {
    }

    int numParameters() const override {
        return vertices.size();
    }

    bool evaluate(double const* params,
                  double* cost,
                  double* jacobians) const override{
        *cost = 0.;
        for (int i = 0; i < vertices.size(); ++i){
            *cost += f.eval(params[i]) * integralOperator[i];
            if(jacobians) jacobians[i] = f.grad(params[i]) * integralOperator[i];
        }

        return true;
    }

    F f;
    Cr::Containers::ArrayView<const Mg::Vector3ui> triangles;
    Cr::Containers::ArrayView<const Mg::Vector3d> vertices;
    Cr::Containers::Array<Mg::Double> integralOperator;
};


template<class F, class L, FunctionalType type>
struct AreaRegularizer : IntegralFunctional<F, L, type> {

    AreaRegularizer(
            Cr::Containers::ArrayView<const Mg::Vector3d> const& vs,
            Cr::Containers::ArrayView<const Mg::Vector3ui> const& ts) :
        IntegralFunctional<F, L, type>(vs, ts)
    {
        auto areas = computeAreas(ts, vs);
        area = std::accumulate(areas.begin(), areas.end(), 0.);
    }

    bool evaluate(double const* params,
                  double* cost,
                  double* jacobians) const override{
        IntegralFunctional<F, L, type>::evaluate(params, &currentArea, jacobians);
        *cost = currentArea - areaRatio * area;
        return true;
    }

    mutable Mg::Double currentArea;
    Mg::Double areaRatio = 0.5, area;
};

using AreaRegularizer1 = AreaRegularizer<Indicator<Mg::Double>, QuadrationLoss, FunctionalType::Area1>;
using AreaRegularizer2 = AreaRegularizer<SmoothStep<Mg::Double>, QuadrationLoss, FunctionalType::Area2>;
using DoubleWellPotential = IntegralFunctional<DoubleWell<Mg::Double>, TrivialLoss, FunctionalType::DoubleWellPotential>;
