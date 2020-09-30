//
// Created by janos on 8/25/20.
//

#include "WeakYamabe.h"
#include "Mesh.h"
#include "C1Functions.h"
#include "Functional.hpp"

#include <Corrade/Containers/Array.h>

#include <adolc/adouble.h>
#include <Eigen/SparseCore>
#include <Eigen/UmfPackSupport>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseQR>
#include <Magnum/Math/Matrix3.h>

#include <igl/gaussian_curvature.h>
#include <igl/cotmatrix.h>

namespace Phasefield {

WeakYamabe::WeakYamabe(Mesh& m) : mesh(&m) {
    mesh->requireGaussianCurvature();
    mesh->requireIntegralOperator();
}

/**
 * ∫ ∇s(∇χ(u)⋅ϕ + χ(u)⋅∇ϕ) - K⋅χ(u)⋅ϕ
 */
template<class Scalar>
void WeakYamabe::operator()(
         ArrayView<const Scalar> parameters, ArrayView<const Scalar> weights, Scalar& out,
         ArrayView<Scalar> gradP, ArrayView<Scalar> gradW) {

    SmoothIndicatorFunction chi;

    size_t n = mesh->vertexCount();
    auto phasefield = parameters.prefix(n);
    auto scalingFactor = parameters.slice(n, 2*n);

    VertexData<Scalar> difference{n};

    for(Face face : mesh->faces()) {
        Math::Vector3<Scalar> gradScaling{0};
        Math::Vector3<Scalar> gradChiOfU{0};

        for(HalfEdge he : face.halfEdges()) {
            Vertex v = he.next().tip();
            Math::Vector3<Scalar> gradBasis{mesh->gradient[he]};
            gradScaling += scalingFactor[v.idx]*gradBasis;
            gradChiOfU += chi.eval(phasefield[v.idx])*gradBasis;
        }

        Scalar scaleDotChi = Math::dot(gradScaling, gradChiOfU);

        for(HalfEdge he : face.halfEdges()) {
            Vertex v = he.next().tip();
            size_t idx = v.idx;

            Scalar chiOfU = chi.eval(phasefield[v.idx]);
            Math::Vector3<Scalar> gradBasis{mesh->gradient[he]};

            difference[v] += (scaleDotChi + chiOfU*Math::dot(gradScaling, gradBasis) - v.gaussianCurvature()*chiOfU)*1./3.*face.area();
        }
    }

    for(Scalar v : difference) {
        out += v*v;
    }
}

size_t WeakYamabe::numParameters() const { return mesh->vertexCount(); }

DEFINE_FUNCTIONAL_CONSTRUCTOR(WeakYamabe)
DEFINE_FUNCTIONAL_OPERATOR(WeakYamabe, double)
DEFINE_FUNCTIONAL_OPERATOR(WeakYamabe, adouble)

struct SmootherStepParametric {
    double edge0, edge1;

    template<class T>
    [[nodiscard]] constexpr T eval(T x) const {
        // Scale, and clamp x to 0..1 range
        x = Math::clamp((x - edge0) / (edge1 - edge0), 0.0, 1.0);
        // Evaluate polynomial
        return x * x * x * (x * (x * 6. - 15.) + 10.);
    }
};


void solveDiffuseYamabeEquation(Mesh& mesh, VertexDataView<const double> parameters, VertexDataView<const double> weights, VertexData<double>& solution) {

    size_t n = mesh.vertexCount();
    size_t m = mesh.faceCount();


    //Eigen::MatrixXi F(m, 3);
    //Eigen::MatrixXd V(n, 3);

    //for(size_t i = 0; i < n; ++i) {
    //    Vector3 p = mesh.positions()[i];
    //    for(size_t j = 0; j < 3; ++j) {
    //        V(i, j) = p[j];
    //    }
    //}

    //size_t i = 0;
    //for(Face face : mesh.faces()) {
    //    size_t j = 0;
    //    for(Vertex v : face.vertices()) {
    //        long idx = v.idx;
    //        assert(idx < n && 0 <= idx);
    //        F(i, j++)  = v.idx;
    //    }
    //    i++;
    //}

    constexpr SmootherStepParametric chi{-1, 1};

    static_assert(chi.eval(-1.) == 0.);
    static_assert(chi.eval(1.) == 1.);

    mesh.requireIntegralOperator();
    mesh.requireGradientOperator();
    mesh.requireGaussianCurvature();
    mesh.requireFaceInformation();

    Array<Eigen::Triplet<double>> triplets;

    for(Face face : mesh.faces()) {

        for(HalfEdge he1 : face.halfEdges()) {
            Vertex v = he1.next().tip();
            Vector3d gradChi = mesh.gradient[he1];
            for(HalfEdge he2 : face.halfEdges()) {
                Vertex w = he2.next().tip();
                double value = Math::dot(mesh.gradient[he2]*chi.eval(parameters[w]), gradChi)*face.area();
                arrayAppend(triplets, InPlaceInit, w.idx, v.idx, value);
            }
        }

        for(Vertex v : face.vertices()) {
            double value = face.area()*Math::pow<2>(1. - chi.eval(parameters[v]));
            arrayAppend(triplets, InPlaceInit, v.idx, v.idx, value);
        }
    }

    Eigen::SparseMatrix<double> A(n, n);
    A.setFromTriplets(triplets.begin(), triplets.end());

    //Eigen::SparseMatrix<double> L;
    //igl::cotmatrix(V, F, L);

    Eigen::UmfPackLU<Eigen::SparseMatrix<double>> solver;
    solver.compute(A);
    if(solver.info() != Eigen::Success) {
        Debug{} << "Umfpack failed";
        switch(solver.info()) {
            case Eigen::NumericalIssue:
                Debug{} << "Numerical Issue";
                break;
            case Eigen::NoConvergence:
                Debug{} << "No Convergence";
                break;
            case Eigen::InvalidInput:
                Debug{} << "Invalid Input";
                break;
        }
    }

    Eigen::VectorXd b = Eigen::VectorXd::Zero(n);
    for(Vertex v : mesh.vertices()) {
        b[v.idx] = mesh.gaussianCurvature[v]*mesh.integral[v]*chi.eval(parameters[v]);
    }

    //Eigen::VectorXd K(n);
    //igl::gaussian_curvature(V, F, K);

    arrayResize(solution, n);
    Eigen::Map<Eigen::VectorXd> map(solution.data(), n);

    map = solver.solve(b);
}

}
