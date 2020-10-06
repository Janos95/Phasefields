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

constexpr Quadratic chi;

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


    static_assert(chi.eval(-1.) == 1.);
    static_assert(chi.eval(1.) == 1.);

    std::vector<Eigen::Triplet<Scalar>> triplets;
    triplets.reserve(12*m);

    {
        for(Face face : mesh.faces()) {

            for(HalfEdge he1 : face.halfEdges()) {
                Vertex v = he1.next().tip();
                if(v.onBoundary()) continue;
                Vec3 gradChi = chi.eval(parameters[v.idx])*Vec3{mesh.gradient[he1]};
                for(HalfEdge he2 : face.halfEdges()) {
                    Vertex w = he2.next().tip();
                    if(w.onBoundary()) continue;
                    Scalar value = Math::dot(Vec3{mesh.gradient[he2]}, gradChi)*face.area();
                    triplets.emplace_back(v.idx, w.idx, value);
                }
            }

            for(Vertex v : face.vertices()) {
                if(!v.onBoundary()) {
                    Scalar value = lambda*face.area()/3.*(1. - chi.eval(parameters[v.idx]));
                    triplets.emplace_back(v.idx, v.idx, value);
                }
            }
        }
    }

    b.resize(n);
    for(Vertex v : mesh.vertices()) {
        if(!v.onBoundary())
            b[v.idx] = 100*mesh.gaussianCurvature[v]*chi.eval(parameters[v])*mesh.integral[v];
        else {
            b[v.idx] = 0;
            triplets.emplace_back(v.idx, v.idx, 1);
        }
    }

    A.resize(n,n);
    A.setFromTriplets(triplets.begin(), triplets.end());

    VecX chiOfU{n}, negChiOfU{n};
    for(size_t j = 0; j < n; ++j) {
        chiOfU[j] = chi.eval(parameters[j]);
        negChiOfU[j] = 1 - chi.eval(parameters[j]);
    }

    //A = chiOfU.asDiagonal()*mesh.fem->stiffness.cast<Scalar>() + lambda*negChiOfU.asDiagonal()*mesh.fem->mass.cast<Scalar>();

    //if constexpr (std::is_same_v<Scalar, double>)
    //    Debug{} << (Atest - A).norm();

    if constexpr(!std::is_same_v<Scalar, adouble>) {
        Eigen::UmfPackLU<Eigen::SparseMatrix<Scalar>> solver{A};
        handleSolverInfo(solver.info());
        x = solver.solve(b);
    } else {
        Eigen::SparseLU<Eigen::SparseMatrix<Scalar>> solver{A};
        handleSolverInfo(solver.info());
        x = solver.solve(b);
    }

}

template<class Scalar>
void assembleAndSolve2(
        Mesh& mesh,
        Scalar lambda,
        VertexDataView<const Scalar> parameters,
        VertexDataView<const Scalar> weights,
        Eigen::SparseMatrix<Scalar>& A1,
        Eigen::SparseMatrix<Scalar>& A2,
        Eigen::Matrix<Scalar, Eigen::Dynamic, 1>& x1,
        Eigen::Matrix<Scalar, Eigen::Dynamic, 1>& x2,
        Eigen::Matrix<Scalar, Eigen::Dynamic, 1>& b1,
        Eigen::Matrix<Scalar, Eigen::Dynamic, 1>& b2) {


    constexpr QuadraticChi chi1;

    size_t n = mesh.vertexCount();
    size_t m = mesh.faceCount();
    using Vec3 = Math::Vector3<Scalar>;
    using VecX = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;

    std::vector<Eigen::Triplet<Scalar>> triplets1;
    std::vector<Eigen::Triplet<Scalar>> triplets2;

    triplets1.reserve(12*m);
    triplets2.reserve(12*m);

    {
        for(Face face : mesh.faces()) {

            for(HalfEdge he1 : face.halfEdges()) {
                Vertex v = he1.next().tip();
                if(v.onBoundary()) continue;
                Vec3 gradChi1 = chi1.eval(parameters[v.idx])*Vec3{mesh.gradient[he1]};
                Vec3 gradChi2 = (1. - chi1.eval(parameters[v.idx]))*Vec3{mesh.gradient[he1]};
                for(HalfEdge he2 : face.halfEdges()) {
                    Vertex w = he2.next().tip();
                    if(w.onBoundary()) continue;
                    Scalar value1 = Math::dot(Vec3{mesh.gradient[he2]}, gradChi1)*face.area();
                    Scalar value2 = Math::dot(Vec3{mesh.gradient[he2]}, gradChi2)*face.area();
                    triplets1.emplace_back(v.idx, w.idx, value1);
                    triplets2.emplace_back(v.idx, w.idx, value2);
                }
            }

            for(Vertex v : face.vertices()) {
                if(!v.onBoundary()) {
                    Scalar value1 = lambda*face.area()/3.*(1. - chi1.eval(parameters[v.idx]));
                    Scalar value2 = lambda*face.area()/3.*chi1.eval(parameters[v.idx]);
                    triplets1.emplace_back(v.idx, v.idx, value1);
                    triplets2.emplace_back(v.idx, v.idx, value2);
                }
            }
        }
    }

    b1.resize(n);
    b2.resize(n);

    for(Vertex v : mesh.vertices()) {
        if(!v.onBoundary()) {
            b1[v.idx] = mesh.gaussianCurvature[v]*chi1.eval(parameters[v])*mesh.integral[v];
            b2[v.idx] = mesh.gaussianCurvature[v]*(1-chi1.eval(parameters[v]))*mesh.integral[v];
        }
        else {
            b1[v.idx] = 0;
            b2[v.idx] = 0;
            triplets1.emplace_back(v.idx, v.idx, 1);
            triplets2.emplace_back(v.idx, v.idx, 1);
        }
    }

    A1.resize(n, n);
    A2.resize(n, n);

    A1.setFromTriplets(triplets1.begin(), triplets1.end());
    A2.setFromTriplets(triplets2.begin(), triplets2.end());

    VecX chiOfU{n}, negChiOfU{n};
    for(size_t j = 0; j < n; ++j) {
        chiOfU[j] = chi.eval(parameters[j]);
        negChiOfU[j] = 1 - chi.eval(parameters[j]);
    }

    //A = chiOfU.asDiagonal()*mesh.fem->stiffness.cast<Scalar>() + lambda*negChiOfU.asDiagonal()*mesh.fem->mass.cast<Scalar>();

    //if constexpr (std::is_same_v<Scalar, double>)
    //    Debug{} << (Atest - A).norm();

    if constexpr(!std::is_same_v<Scalar, adouble>) {
        Eigen::UmfPackLU<Eigen::SparseMatrix<Scalar>> solver;
        solver.solve(A1);
        handleSolverInfo(solver.info());
        x1 = solver.solve(b1);

        solver.solve(A2);
        handleSolverInfo(solver.info());
        x2 = solver.solve(b2);
    } else {
        Eigen::SparseLU<Eigen::SparseMatrix<Scalar>> solver;
        solver.solve(A1);
        handleSolverInfo(solver.info());
        x1 = solver.solve(b1);

        solver.solve(A2);
        handleSolverInfo(solver.info());
        x2 = solver.solve(b2);
    }
}

    template<class Scalar>
void DiffuseYamabe::operator()(
         ArrayView<const Scalar> parameters, ArrayView<const Scalar> weights, Scalar& out,
         ArrayView<Scalar> gradP, ArrayView<Scalar> gradW) {

    using Vec3 = Math::Vector3<Scalar>;
    using VecX = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;

    size_t n = mesh.vertexCount();
    size_t m = mesh.faceCount();
    double lambda = lambdaWeight*(lambdaScaling ? *lambdaScaling : 1.);

    Eigen::SparseMatrix<Scalar> A;
    VecX b,x;
    assembleAndSolve<Scalar>(mesh, lambda, parameters, weights, A, x, b);

    VecX gradX{gradP ? n : 0};
    gradX.setZero();

    if(energy == EnergyType::Dirichlet) {
        for(Face face : mesh.faces()) {
            Math::Vector3<Scalar> grad{0};
            Scalar rescaling = (curvatureRescaling ? getRescalingFactor(face) : 1.);

            for(HalfEdge he : face.halfEdges()) {
                Vertex v = he.next().tip();
                grad += x[v.idx]*Vec3{mesh.gradient[he]};
            }

            if(gradP) {
                for(HalfEdge he : face.halfEdges()){
                    Vertex v = he.next().tip();
                    gradX[v.idx] += rescaling*mesh.faceArea[face]*2*Math::dot(grad, Vec3{mesh.gradient[he]});
                }
            }

            Scalar gradNormSquared = grad.dot() * rescaling;
            //CORRADE_INTERNAL_ASSERT(!Math::isNan(gradNormSquared));
            out += mesh.faceArea[face]*gradNormSquared;
        }
    } else if(energy == EnergyType::Hencky) {
        for(Vertex v : mesh.vertices()) {
            Scalar rescaling = (curvatureRescaling ? getRescalingFactor(v) : 1.);
            out += x[v.idx]*x[v.idx]*mesh.integral[v];
            if(gradP) {
                gradX[v.idx] += 2*x[v.idx]*mesh.integral[v]*rescaling;
                //gradP[v.idx] += x[v.idx]*x[v.idx]*mesh.integral[v]*chi.grad(parameters[v.idx]);
            }
        }
    }

    if(gradP) {
        if constexpr (!std::is_same_v<Scalar, adouble>) {

            Eigen::UmfPackLU<Eigen::SparseMatrix<Scalar>> solver{A.transpose()};
            VecX l = solver.solve(gradX);

            VecX bp{n};
            for(Vertex v : mesh.vertices()) {
                if(!v.onBoundary())
                    bp[v.idx] = 100*mesh.gaussianCurvature[v]*mesh.integral[v]*chi.grad(parameters[v.idx]);
                else
                    bp[v.idx] = 0;
            }

            std::vector<Eigen::Triplet<Scalar>> triplets;
            triplets.reserve(12*m);

            {
                for(Face face : mesh.faces()) {

                    for(HalfEdge he1 : face.halfEdges()) {
                        Vertex v = he1.next().tip();
                        if(v.onBoundary()) continue;
                        auto gradChi = chi.grad(parameters[v.idx])*Vec3{mesh.gradient[he1]};
                        for(HalfEdge he2 : face.halfEdges()) {
                            Vertex w = he2.next().tip();
                            if(w.onBoundary()) continue;
                            Scalar value = Math::dot(Vec3{mesh.gradient[he2]}, gradChi)*face.area();
                            triplets.emplace_back(v.idx, w.idx, value);
                        }
                    }

                    for(Vertex v : face.vertices()) {
                        if(!v.onBoundary()) {
                            Scalar value = lambda*face.area()/3.*(-chi.grad(parameters[v.idx]));
                            triplets.emplace_back(v.idx, v.idx, value);
                        }
                    }
                }
            }

            Eigen::SparseMatrix<Scalar> Ap(n,n);
            Ap.setFromTriplets(triplets.begin(), triplets.end());

            //Eigen::SparseMatrix<Scalar> Ap = chiPrime.asDiagonal()*mesh.fem->stiffness - lambda*chiPrime.asDiagonal()*mesh.fem->mass;
            Eigen::Map<VecX>{gradP.data(), long(gradP.size())} -= l.asDiagonal()*(Ap*x - bp);
            Debug{} << "Gradient norm : " << (l.asDiagonal()*(Ap*x-bp)).norm();
            Debug{} << "l norm : " << l.norm();
            Debug{} << "Ap*x - bp norm : " << (Ap*x - bp).norm();
            Debug{} << "Ap - Ap norm : " << (Ap - A).norm();
            Debug{} << "bp - b norm : " << (bp - b).norm();
            Debug{} << "grad x norm" << gradX.norm();
            Debug{} << "x norm" << x.norm();
        }
    }
}

size_t DiffuseYamabe::numParameters() const { return mesh.vertexCount(); }

void DiffuseYamabe::drawImGuiOptions(VisualizationProxy& proxy) {
    bool redraw = false;
    if(ImGui::Checkbox("Draw Solution", &drawSolution)) {
        if(drawSolution) {
            proxy.setCallbacks(
                    [this, &proxy](Node node) {
                        Eigen::SparseMatrix<double> A;
                        Eigen::VectorXd b,x;
                        double lambda = lambdaWeight*(lambdaScaling ? *lambdaScaling : 1.);
                        assembleAndSolve<double>(mesh, lambda, node.phasefield(), node.temporary(), A, x, b);
                        ArrayView<double> solution {x.data(), size_t(x.size())};
                        proxy.drawValuesNormalized(solution);
                    },
                    [this]{ drawSolution = false; });
        } else proxy.setDefaultCallback();
        redraw = true;
    }

    ImGui::SameLine();

    if(ImGui::Checkbox("Draw Solution Thresholded", &drawSolutionThresholded)) {
        if(drawSolutionThresholded) {
            proxy.setCallbacks(
                    [this, &proxy](Node node) {
                        auto phasefield = node.phasefield();
                        size_t n = phasefield.size();

                        VertexData<bool> inInterior{DirectInit, n, 1};
                        for(Edge e : mesh.edges()) {
                            Vertex v1 = e.vertex1();
                            Vertex v2 = e.vertex2();
                            if(Math::sign(phasefield[v1]) != Math::sign(phasefield[v2])) {
                                inInterior[v1] = false;
                                inInterior[v2] = false;
                            }
                            for(Vertex v : {v1, v2}) {
                                if(inInterior[v] && v.onBoundary())
                                    inInterior[v1] = false;
                            }
                        }

                        Array<Eigen::Triplet<double>> triplets;
                        Eigen::SparseMatrix<double> A(n,n);
                        Eigen::VectorXd b(n);

                        for(Vertex v : mesh.vertices()) {
                            if(inInterior[v])
                                b[v.idx] = mesh.gaussianCurvature[v]*mesh.integral[v]; //not really true of the boundary ..
                            else {
                                b[v.idx] = 0;
                                arrayAppend(triplets, InPlaceInit, v.idx, v.idx, 1.);
                            }
                        }


                        for(Face face : mesh.faces()) {
                            for(HalfEdge he1 : face.halfEdges()) {
                                Vertex v = he1.next().tip();
                                if(!inInterior[v]) continue;
                                for(HalfEdge he2 : face.halfEdges()) {
                                    Vertex w = he2.next().tip();
                                    if(!inInterior[v]) continue;
                                    double value = Math::dot(mesh.gradient[he2], mesh.gradient[he1])*face.area();
                                    arrayAppend(triplets, InPlaceInit, v.idx, w.idx, value);
                                }
                            }
                        }

                        A.setFromTriplets(triplets.begin(), triplets.end());

                        Eigen::UmfPackLU<Eigen::SparseMatrix<double>> solver{A};
                        handleSolverInfo(solver.info());
                        Eigen::VectorXd x = solver.solve(b);

                        Debug{} << "Solution norm" << x.norm();

                        proxy.drawValuesNormalized({x.data(), size_t(x.size())});
                    },
                    [this]{ drawSolutionThresholded = false; });
        } else proxy.setDefaultCallback();
        redraw = true;
    }

    if(ImGui::Checkbox("Draw Solution Gradient", &drawSolutionGradient)) {
        if(drawSolutionGradient) {
            proxy.setCallbacks(
                    [this, &proxy](Node node) {
                        Eigen::SparseMatrix<double> A;
                        Eigen::VectorXd b,x;
                        double lambda = lambdaWeight*(lambdaScaling ? *lambdaScaling : 1.);
                        assembleAndSolve<double>(mesh, lambda , node.phasefield(), node.temporary(), A, x, b);

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
        redraw = true;
    }

    static const double minWeight = 0.0001;
    static const double maxWeight = 100;
    if(ImGui::DragScalar("lambda", ImGuiDataType_Double, &lambdaWeight, 1.f, &minWeight, &maxWeight, "%f", 1)) {
        redraw |= drawSolutionThresholded || drawSolution || drawSolutionGradient;
    }


    if(ImGui::BeginCombo("##energie", EnergyType::to_string(energy))) {
        for(auto type : EnergyType::range) {
            bool isSelected = type == energy;
            if(ImGui::Selectable(EnergyType::to_string(type), isSelected))
                energy = type;
            if(isSelected)
                ImGui::SetItemDefaultFocus();
        }
        ImGui::EndCombo();
    }

    if(ImGui::Checkbox("Curvature Rescaling", &curvatureRescaling)) {
        redraw |= drawSolutionThresholded || drawSolution || drawSolutionGradient;
    }

    if(redraw) proxy.redraw();
}

double DiffuseYamabe::getRescalingFactor(Face f) const {
    double c = 0;
    constexpr double eps = 1e-7;
    for(Vertex v : f.vertices()) {
        c += mesh.gaussianCurvature[v];
    }
    c /= 3.;
    double r = 1./(sqrt(c) + eps);
    return r*r*r;
}

double DiffuseYamabe::getRescalingFactor(Vertex v) const {
    double c = 0;
    constexpr double eps = 1e-7;
    double r = 1./(sqrt(mesh.gaussianCurvature[v]) + eps);
    return r*r*r;
}

DEFINE_FUNCTIONAL_CONSTRUCTOR(DiffuseYamabe)
DEFINE_FUNCTIONAL_OPERATOR(DiffuseYamabe, double)
DEFINE_FUNCTIONAL_OPERATOR(DiffuseYamabe, adouble)

}