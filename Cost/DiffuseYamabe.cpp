//
// Created by janos on 8/25/20.
//

#include "DiffuseYamabe.h"
#include "Mesh.h"
#include "C1Functions.h"
#include "Tree.h"
#include "Functional.hpp"
#include "FEM.hpp"
#include "VisualizationProxy.h"
#include <ScopedTimer/ScopedTimer.h>

#include <Corrade/Containers/Array.h>

#include <adolc/adouble.h>
#include <Eigen/SparseCore>
#include <Eigen/UmfPackSupport>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseLU>
//#include <eigen3/unsupported/Eigen/AdolcForward>
#include <Eigen/SparseCholesky>
#include <Magnum/Math/Matrix3.h>

#include <imgui.h>

namespace Eigen {

template<> struct NumTraits<adouble>
        : NumTraits<double> // permits to get the epsilon, dummy_precision, lowest, highest functions
{
    typedef adouble Real;
    typedef adouble NonInteger;
    typedef adouble Nested;

    enum {
        IsComplex = 0,
        IsInteger = 0,
        IsSigned = 1,
        RequireInitialization = 1,
        ReadCost = 1,
        AddCost = 3,
        MulCost = 3
    };
};

template<> struct NumTraits<adub>
        : NumTraits<double> // permits to get the epsilon, dummy_precision, lowest, highest functions
{
    typedef adouble Real;
    typedef adouble NonInteger;
    typedef adouble Nested;

    enum {
        IsComplex = 0,
        IsInteger = 0,
        IsSigned = 1,
        RequireInitialization = 1,
        ReadCost = 1,
        AddCost = 3,
        MulCost = 3
    };
};

}

inline const adouble& conj(const adouble& x)  { return x; }
inline const adouble& real(const adouble& x)  { return x; }
inline adouble imag(const adouble&)    { return 0.; }
inline adouble abs(const adouble&  x)  { return fabs(x); }
inline adouble abs2(const adouble& x)  { return x*x; }

namespace Phasefield {

#define CHI QuadraticChi

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

void handleSolverInfo(Eigen::ComputationInfo info) {
    switch(info) {
        case Eigen::NumericalIssue:
            Debug{} << "Numerical Issue";
            break;
        case Eigen::NoConvergence:
            Debug{} << "No Convergence";
            break;
        case Eigen::InvalidInput:
            Debug{} << "Invalid Input";
            break;
        case Eigen::Success:
            Debug{} << "Solver Successfull";
            break;
    }
}

DiffuseYamabe::DiffuseYamabe(Mesh& m) : mesh(m) {
    mesh.requireIntegralOperator();
    mesh.requireGradientOperator();
    mesh.requireGaussianCurvature();
    mesh.requireFaceInformation();
    mesh.requireStiffnessMatrix();
    mesh.requireMassMatrix();
}

template<class Scalar>
void assembleAndSolve(
        Mesh& mesh,
        Scalar lambda,
        VertexDataView<const Scalar> parameters,
        VertexDataView<const Scalar> weights,
        Eigen::SparseMatrix<Scalar>& A,
        Eigen::Matrix<Scalar, Eigen::Dynamic, 1>& x,
        Eigen::Matrix<Scalar, Eigen::Dynamic, 1>& b) {

    size_t n = mesh.vertexCount();
    size_t m = mesh.faceCount();
    using Vec3 = Math::Vector3<Scalar>;
    using VecX = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;

    constexpr CHI chi;

    static_assert(chi.eval(-1.) == 0.);
    static_assert(chi.eval(1.) == 1.);

    mesh.requireIntegralOperator();
    mesh.requireGradientOperator();
    mesh.requireGaussianCurvature();
    mesh.requireFaceInformation();
    mesh.requireStiffnessMatrix();
    mesh.requireMassMatrix();

    //Array<Eigen::Triplet<Scalar>> triplets{12*m};

    //{
    //    size_t i = 0;
    //    for(Face face : mesh.faces()) {

    //        for(HalfEdge he1 : face.halfEdges()) {
    //            Vertex v = he1.next().tip();
    //            auto gradChi = chi.eval(parameters[v.idx])*Vec3{mesh.gradient[he1]};
    //            for(HalfEdge he2 : face.halfEdges()) {
    //                Vertex w = he2.next().tip();
    //                Scalar value = Math::dot(Vec3{mesh.gradient[he2]}, gradChi)*face.area();
    //                triplets[i++] = Eigen::Triplet<Scalar>{v.idx, w.idx, value};
    //            }
    //        }

    //        for(Vertex v : face.vertices()) {
    //            Scalar value = face.area()/3.*(1. - chi.eval(parameters[v.idx]));
    //            CORRADE_ASSERT(i < triplets.size(), "Index out of bounds",);
    //            triplets[i++] = Eigen::Triplet<Scalar>{v.idx, v.idx, value};
    //        }

    //    }
    //}

    //A.resize(n,n);
    //A.setFromTriplets(triplets.begin(), triplets.end());

    VecX chiOfU{n}, negChiOfU{n};
    for(size_t j = 0; j < n; ++j) {
        chiOfU[j] = chi.eval(parameters[j]);
        negChiOfU[j] = 1 - chi.eval(parameters[j]);
    }

    A = chiOfU.asDiagonal()*mesh.fem->stiffness.cast<Scalar>() + lambda*negChiOfU.asDiagonal()*mesh.fem->mass.cast<Scalar>();

    //if constexpr (std::is_same_v<Scalar, double>)
    //    Debug{} << (Atest - A).norm();

    b.resize(n);
    for(Vertex v : mesh.vertices()) {
        b[v.idx] = mesh.gaussianCurvature[v]*mesh.integral[v]*chi.eval(parameters[v.idx]);
    }

    if constexpr(!std::is_same_v<Scalar, adouble>) {
        ScopedTimer t{"Umfpack", true};
        Eigen::UmfPackLU<Eigen::SparseMatrix<Scalar>> solver{A};
        handleSolverInfo(solver.info());
        x = solver.solve(b);
    } else {
        ScopedTimer t{"Eigen LU", true};
        Eigen::SparseLU<Eigen::SparseMatrix<Scalar>> solver{A};
        handleSolverInfo(solver.info());
        x = solver.solve(b);
    }

}

template<class Scalar>
void DiffuseYamabe::operator()(
         ArrayView<const Scalar> parameters, ArrayView<const Scalar> weights, Scalar& out,
         ArrayView<Scalar> gradP, ArrayView<Scalar> gradW) {

    ScopedTimer t{"Yamabe", true};

    adouble s,r;
    s *= r;

    using Vec3 = Math::Vector3<Scalar>;
    using VecX = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;

    size_t n = mesh.vertexCount();
    size_t m = mesh.faceCount();

    Eigen::SparseMatrix<Scalar> A;
    VecX b,x;
    assembleAndSolve<Scalar>(mesh, lambda, parameters, weights, A, x, b);

    VecX gradDirichlet{gradP ? n : 0};
    gradDirichlet.setZero();

    for(Face face : mesh.faces()) {
        Math::Vector3<Scalar> grad{0};

        for(HalfEdge he : face.halfEdges()) {
            Vertex v = he.next().tip();
            grad += x[v.idx]*Vec3{mesh.gradient[he]};
        }

        if(gradP) {
            for(HalfEdge he : face.halfEdges()){
                Vertex v = he.next().tip();
                gradDirichlet[v.idx] += mesh.faceArea[face]*2*Math::dot(grad, Vec3{mesh.gradient[he]});
            }
        }

        Scalar gradNormSquared = grad.dot();
        //CORRADE_INTERNAL_ASSERT(!Math::isNan(gradNormSquared));
        out += mesh.faceArea[face]*gradNormSquared;
    }

    if(gradP) {
        if constexpr (!std::is_same_v<Scalar, adouble>) {
            constexpr CHI chi;

            Eigen::UmfPackLU<Eigen::SparseMatrix<Scalar>> solver{A.transpose()};
            VecX l = solver.solve(gradDirichlet);

            VecX chiPrime{n}, bp{n};
            for(size_t j = 0; j < n; ++j) {
                Scalar chiGrad = chi.grad(parameters[j]);
                bp[j] = mesh.gaussianCurvature[j]*mesh.integral[j]*chiGrad;
                chiPrime[j] = chiGrad;
            }

            Eigen::SparseMatrix<Scalar> Ap = chiPrime.asDiagonal()*mesh.fem->stiffness - lambda*chiPrime.asDiagonal()*mesh.fem->mass;
            Eigen::Map<VecX>{gradP.data(), long(gradP.size())} -= l.asDiagonal()*(Ap*x - bp);
        }
    }
}

size_t DiffuseYamabe::numParameters() const { return mesh.vertexCount(); }

void DiffuseYamabe::drawImGuiOptions(VisualizationProxy& proxy) {
    if(ImGui::Checkbox("Draw Solution", &drawSolution)) {
        if(drawSolution) {
            proxy.setCallbacks(
                    [this, &proxy](Node node) {
                        Eigen::SparseMatrix<double> A;
                        Eigen::VectorXd b,x;
                        assembleAndSolve<double>(mesh, lambda, node.phasefield(), node.temporary(), A, x, b);
                        ArrayView<double> solution {x.data(), size_t(x.size())};
                        proxy.drawValuesNormalized(solution);
                    },
                    [this]{ drawSolution = false; });
        } else proxy.setDefaultCallback();
        proxy.redraw();
    }

    ImGui::SameLine();

    if(ImGui::Checkbox("Draw Solution Thresholded", &drawSolutionThresholded)) {
        if(drawSolutionThresholded) {
            proxy.setCallbacks(
                    [this, &proxy](Node node) {
                        auto phasefield = node.phasefield();
                        Array<Face> facesInPhase;
                        Array<Vertex> verticesInPhase;
                        VertexData<size_t> vertexIndices{DirectInit, phasefield.size(), Invalid};

                        for(Face f : mesh.faces()) {
                            double average = 0;
                            for(Vertex v : f.vertices()) {
                                average += phasefield[v];
                            }
                            if(average > 0) {
                                arrayAppend(facesInPhase, f);
                                for(Vertex v : f.vertices()) {
                                    if(vertexIndices[v] == Invalid) {
                                        arrayAppend(verticesInPhase, v);
                                        vertexIndices[v] = verticesInPhase.size() - 1;
                                    }
                                }
                            }
                        }

                        size_t n = verticesInPhase.size();

                        Eigen::SparseMatrix<double> A(n,n);
                        Eigen::VectorXd b(n);

                        for(size_t i = 0; i < n; ++i) {
                            Vertex v = verticesInPhase[i];
                            b[i] = mesh.gaussianCurvature[v]*mesh.integral[v]; //not really true of the boundary ..
                        }

                        Array<Eigen::Triplet<double>> triplets;

                        for(Face face : facesInPhase) {
                            for(HalfEdge he1 : face.halfEdges()) {
                                Vertex v = he1.next().tip();
                                size_t vIdx = vertexIndices[v];
                                for(HalfEdge he2 : face.halfEdges()) {
                                    Vertex w = he2.next().tip();
                                    double value = Math::dot(mesh.gradient[he2], mesh.gradient[he1])*face.area();
                                    size_t wIdx = vertexIndices[w];
                                    arrayAppend(triplets, InPlaceInit, vIdx, wIdx, value);
                                }
                            }
                        }

                        A.setFromTriplets(triplets.begin(), triplets.end());

                        Eigen::UmfPackLU<Eigen::SparseMatrix<double>> solver{A};
                        handleSolverInfo(solver.info());
                        Eigen::VectorXd x = solver.solve(b);

                        double min = x.minCoeff();
                        double max = x.maxCoeff();

                        for(size_t i = 0; i < n; ++i) {
                            x[i] = (x[i] - min)/(max - min);
                        }

                        VertexData<double> solution{phasefield.size()};
                        for(Vertex v : verticesInPhase) {
                            size_t i = vertexIndices[v];
                            solution[v] = x[i];
                        }
                        proxy.drawValues(solution);
                    },
                    [this]{ drawSolutionThresholded = false; });
        } else proxy.setDefaultCallback();
        proxy.redraw();
    }

    if(ImGui::Checkbox("Draw Solution Gradient", &drawSolutionGradient)) {
        if(drawSolutionGradient) {
            proxy.setCallbacks(
                    [this, &proxy](Node node) {
                        Eigen::SparseMatrix<double> A;
                        Eigen::VectorXd b,x;
                        assembleAndSolve<double>(mesh, lambda, node.phasefield(), node.temporary(), A, x, b);

                        VertexData<double> gradSolution{node.phasefield().size()};

                        for(Face face : mesh.faces()) {
                            Vector3d grad;

                            for(HalfEdge he : face.halfEdges()) {
                                Vertex v = he.next().tip();
                                grad += x[v.idx]*mesh.gradient[he];
                            }

                            double gradNorm = grad.length();

                            for(Vertex v : face.vertices()) {
                                gradSolution[v] += gradNorm;
                            }
                        }

                        for(Vertex v : mesh.vertices())
                            gradSolution[v] /= double(mesh.degree[v]);

                        proxy.drawValuesNormalized(gradSolution);
                    },
                    [this]{ drawSolutionGradient = false; });
        } else proxy.setDefaultCallback();
        proxy.redraw();
    }

    static const double minWeight = 0.0001;
    static const double maxWeight = 100;
    ImGui::DragScalar("lambda", ImGuiDataType_Double, &lambda, 1.f, &minWeight, &maxWeight, "%f", 1);
}

DEFINE_FUNCTIONAL_CONSTRUCTOR(DiffuseYamabe)
DEFINE_FUNCTIONAL_OPERATOR(DiffuseYamabe, double)
DEFINE_FUNCTIONAL_OPERATOR(DiffuseYamabe, adouble)

}
