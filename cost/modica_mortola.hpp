//
// Created by janos on 27.11.19.
//

#pragma once

#include "functional.hpp"

#include <Corrade/Containers/ArrayView.h>

namespace Eigen {
    template<typename _Scalar, int _Options, typename _StorageIndex>
    class SparseMatrix;
}

namespace Cr = Corrade;
namespace Mg = Magnum;

struct DirichletEnergy :  Functional
{

    DirichletEnergy(
            Cr::Containers::ArrayView<const Vector3> const& vertices,
            Cr::Containers::ArrayView<const Vector3ui> const& faces,
            Float epsilon);

    int NumParameters() const override;

    bool Evaluate(double const* parameters,
                  double* residual,
                  double* jacobians) const override;

    double m_epsilon;

    Cr::Containers::ArrayView<const Vector3ui> m_triangles;
    Cr::Containers::ArrayView<const Vector3> m_vertices;
    Cr::Containers::Array<Mg::Float> areas;

    Cr::Containers::Pointer<Eigen::SparseMatrix<Mg::Float, 0 /* ColMajor */, int>> m_GSQ;
};


struct DoubleWellPotential :  ceres::FirstOrderFunction
{

    DoubleWellPotential(
        Cr::Containers::ArrayView<const Vector3> const& vertices,
        Cr::Containers::ArrayView<const Vector3ui> const& faces,
        Mg::Float epsilon);

    int NumParameters() const override;

    bool Evaluate(double const* parameters,
                  double* cost,
                  double* jacobians) const override;

    Cr::Containers::ArrayView<const Vector3ui> m_triangles;
    Cr::Containers::ArrayView<const Vector3> m_vertices;
    Cr::Containers::Array<Mg::Float> areas;

    Mg::Float m_epsilon;
};


struct AreaRegularizer : public ceres::FirstOrderFunction
{

    AreaRegularizer(
        Cr::Containers::ArrayView<const Vector3> const& vertices,
        Cr::Containers::ArrayView<const Vector3ui> const& faces,
        Mg::Float areaRatio = .5f);

    int NumParameters() const override;

    bool Evaluate(double const* parameters,
                  double* cost,
                  double* jacobians) const override;


    Cr::Containers::ArrayView<const Vector3ui> triangles;
    Cr::Containers::ArrayView<const Vector3> vertices;
    Cr::Containers::Array<Mg::Float> areas;

    Mg::Float area;
    Mg::Float areaRatio;
    Mg::Float currentArea;
};

