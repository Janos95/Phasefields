//
// Created by janos on 8/25/20.
//

#include "DiffuseYamabe.h"
#include "EigenHelper.h"
#include "Mesh.h"
#include "C1Functions.h"
#include "Tree.h"
#include "Functional.hpp"
#include "FEM.hpp"
#include "VisualizationProxy.h"
#include <ScopedTimer/ScopedTimer.h>

#include <Corrade/Containers/Array.h>
#include <Magnum/Math/Matrix3.h>

#include <imgui.h>

namespace Phasefield {

void getFaceData(Face face, Vertex* vs, HalfEdge* hs) {
    HalfEdge he = face.halfEdge();
    for(size_t i = 0; i < 3; ++i, he = he.next()) {
        hs[i] = he;
        vs[(i+1)%3] = he.tail();
    }
}

DiffuseYamabe::DiffuseYamabe(Mesh& m) : mesh(m) {
    mesh.requireIntegralOperator();
    mesh.requireGradientOperator();
    mesh.requireGaussianCurvature();
    mesh.requireFaceInformation();
    mesh.requireStiffnessMatrix();
    mesh.requireBoundaryInformation();
}

constexpr QuadraticChi chi1;
constexpr QuadraticChiMirrored chi2;

template<class Scalar>
void assembleAndSolve(
        DiffuseYamabe& yamabe,
        Scalar lambda,
        VertexDataView<const Scalar> parameters,
        [[maybe_unused]] VertexDataView<const Scalar> weights,
        Eigen::SparseMatrix<Scalar>& A1,
        SolverType<Scalar>& solver1,
        Eigen::SparseMatrix<Scalar>& A2,
        SolverType<Scalar>& solver2,
        Eigen::Matrix<Scalar, Eigen::Dynamic, 1>& x1,
        Eigen::Matrix<Scalar, Eigen::Dynamic, 1>& x2,
        Eigen::Matrix<Scalar, Eigen::Dynamic, 1>& b1,
        Eigen::Matrix<Scalar, Eigen::Dynamic, 1>& b2) {

    using Vec3 = Math::Vector3<Scalar>;
    using VecX = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;

    Mesh& mesh = yamabe.mesh;

    size_t n = mesh.vertexCount();
    size_t m = mesh.faceCount();

    bool doPositive = yamabe.positivePhase;
    bool doNegative = yamabe.negativePhase;

    std::vector<Eigen::Triplet<Scalar>> triplets1;
    std::vector<Eigen::Triplet<Scalar>> triplets2;

    triplets1.reserve(12*m);
    triplets2.reserve(12*m);

    b1.resize(n);
    b2.resize(n);

    {
        auto& stiffnessElements = mesh.stiffnessElements;

        for(Face face : mesh.faces()) {

            HalfEdge hs[3];
            Vertex vs[3];
            getFaceData(face, vs, hs);

            Scalar uT = 0;

            for(Vertex v : vs) {
                uT += parameters[v.idx];
            }

            uT /= 3.;
            Scalar chi = chi1.eval(uT);

            for(size_t i = 0; i < 3; ++i) {
                Vertex u = vs[i];
                Vertex v = vs[(i+1)%3];
                Vertex w = vs[(i+2)%3];
                if(!v.onBoundary() && !w.onBoundary()) {
                    triplets1.emplace_back(v.idx, w.idx, stiffnessElements[hs[i]]*chi);
                    triplets1.emplace_back(w.idx, v.idx, stiffnessElements[hs[i]]*chi);

                    triplets2.emplace_back(v.idx, w.idx, stiffnessElements[hs[i]]*(1.-chi));
                    triplets2.emplace_back(w.idx, v.idx, stiffnessElements[hs[i]]*(1.-chi));
                }

                if(!u.onBoundary()) {
                    double gradientSquared = hs[i].gradient().dot();
                    triplets1.emplace_back(u.idx, u.idx, gradientSquared*chi*face.area());
                    triplets2.emplace_back(u.idx, u.idx, gradientSquared*(1.-chi)*face.area());
                }
            }
        }

        for(Vertex v : mesh.vertices()) {
            if(!v.onBoundary()) {
                Scalar K = mesh.gaussianCurvature[v];
                Scalar chi = chi1.eval(parameters[v.idx]);

                b1[v.idx] = chi*mesh.integral[v]*K;
                b2[v.idx] = (1-chi)*mesh.integral[v]*K;

                triplets1.emplace_back(v.idx, v.idx, lambda*(1-chi)*mesh.integral[v]);
                triplets2.emplace_back(v.idx, v.idx, lambda*chi*mesh.integral[v]);
            } else {
                b1[v.idx] = 0;
                b2[v.idx] = 0;
                triplets1.emplace_back(v.idx, v.idx, 1);
                triplets2.emplace_back(v.idx, v.idx, 1);
            }
        }
    }

    A1.resize(n, n);
    A2.resize(n, n);


    if(doPositive) {
        A1.setFromTriplets(triplets1.begin(), triplets1.end());
    }
    if(doNegative) {
        A2.setFromTriplets(triplets2.begin(), triplets2.end());
    }

    if(doPositive) {
        solver1.compute(A1);
        handleSolverInfo(solver1.info());
        x1 = solver1.solve(b1);
    }

    if(doNegative) {
        solver2.compute(A2);
        handleSolverInfo(solver2.info());
        x2 = solver2.solve(b2);
    }
}

void shapeDerivative(
        DiffuseYamabe& yamabe,
        VertexDataView<const double> parameters,
        double lambda,
        Eigen::SparseMatrix<double> const& A1,
        SolverType<double>& solver1,
        Eigen::SparseMatrix<double> const& A2,
        SolverType<double>& solver2,
        Eigen::VectorXd const& x1,
        Eigen::VectorXd const& x2,
        Eigen::VectorXd const& gradX1,
        Eigen::VectorXd const& gradX2,
        ArrayView<double> grad){

    Mesh& mesh = yamabe.mesh;
    size_t n = mesh.vertexCount();

    bool doPositive = yamabe.positivePhase;
    bool doNegative = yamabe.negativePhase;

    Eigen::VectorXd l1, l2;

    if(doPositive) {
        l1 = solver1.solve(gradX1);
    }
    if(doNegative) {
        l2 = solver2.solve(gradX2);
    }

    std::vector<Eigen::Triplet<double>> triplets1;
    std::vector<Eigen::Triplet<double>> triplets2;

    triplets1.reserve(12*mesh.faceCount());
    triplets2.reserve(12*mesh.faceCount());

    {
        auto& stiffnessElements = mesh.stiffnessElements;

        for(Face face : mesh.faces()) {

            HalfEdge hs[3];
            Vertex vs[3];
            getFaceData(face, vs, hs);

            double uT = 0;
            double s1T = 0;
            double s2T = 0;

            for(Vertex v : vs) {
                uT += parameters[v.idx];
                if(doPositive)
                    s1T += x1[v.idx];
                if(doNegative)
                    s2T += x2[v.idx];
            }

            uT /= 3.;
            s1T /= 3.;
            s2T /= 3.;

            double chiGrad = chi1.grad(uT)/3.;

            for(size_t i = 0; i < 3; ++i) {
                Vertex u = vs[i];
                Vertex v = vs[(i+1)%3];
                Vertex w = vs[(i+2)%3];
                if(!v.onBoundary() && !w.onBoundary()) {
                    if(doPositive) {
                        triplets1.emplace_back(v.idx, w.idx, stiffnessElements[hs[i]]*chiGrad*x1[w.idx]);
                        triplets1.emplace_back(v.idx, u.idx, stiffnessElements[hs[i]]*chiGrad*x1[w.idx]);
                        triplets1.emplace_back(v.idx, v.idx, stiffnessElements[hs[i]]*chiGrad*x1[w.idx]);

                        triplets1.emplace_back(w.idx, v.idx, stiffnessElements[hs[i]]*chiGrad*x1[v.idx]);
                        triplets1.emplace_back(w.idx, u.idx, stiffnessElements[hs[i]]*chiGrad*x1[v.idx]);
                        triplets1.emplace_back(w.idx, w.idx, stiffnessElements[hs[i]]*chiGrad*x1[v.idx]);
                    }

                    if(doNegative) {
                        triplets2.emplace_back(v.idx, u.idx, stiffnessElements[hs[i]]*(-chiGrad)*x2[w.idx]);
                        triplets2.emplace_back(v.idx, w.idx, stiffnessElements[hs[i]]*(-chiGrad)*x2[w.idx]);
                        triplets2.emplace_back(v.idx, v.idx, stiffnessElements[hs[i]]*(-chiGrad)*x2[w.idx]);

                        triplets2.emplace_back(w.idx, v.idx, stiffnessElements[hs[i]]*(-chiGrad)*x2[v.idx]);
                        triplets2.emplace_back(w.idx, u.idx, stiffnessElements[hs[i]]*(-chiGrad)*x2[v.idx]);
                        triplets2.emplace_back(w.idx, w.idx, stiffnessElements[hs[i]]*(-chiGrad)*x2[v.idx]);
                    }
                }

                if(!u.onBoundary()) {
                    double gradientSquared = hs[i].gradient().dot();
                    if(doPositive) {
                        triplets1.emplace_back(u.idx, u.idx, gradientSquared*chiGrad*face.area()*x1[u.idx]);
                        triplets1.emplace_back(u.idx, v.idx, gradientSquared*chiGrad*face.area()*x1[u.idx]);
                        triplets1.emplace_back(u.idx, w.idx, gradientSquared*chiGrad*face.area()*x1[u.idx]);
                    }
                    if(doNegative) {
                        triplets2.emplace_back(u.idx, u.idx, gradientSquared*(-chiGrad)*face.area()*x2[u.idx]);
                        triplets2.emplace_back(u.idx, v.idx, gradientSquared*(-chiGrad)*face.area()*x2[u.idx]);
                        triplets2.emplace_back(u.idx, w.idx, gradientSquared*(-chiGrad)*face.area()*x2[u.idx]);
                    }
                }
            }
        }

        for(Vertex v : mesh.vertices()) {
            if(!v.onBoundary()) {
                double K = mesh.gaussianCurvature[v];
                double chiGradIntegrated = chi1.grad(parameters[v.idx])*mesh.integral[v];

                if(doPositive) {
                    double vv1 = (-x1[v.idx]*lambda - K)*chiGradIntegrated;
                    triplets1.emplace_back(v.idx, v.idx, vv1);
                }

                if(doNegative) {
                    double vv2 = (x2[v.idx]*lambda + K)*chiGradIntegrated;
                    triplets2.emplace_back(v.idx, v.idx, vv2);
                }
            }
        }
    }

    Eigen::SparseMatrix<double> sensitivity1(n,n);
    Eigen::SparseMatrix<double> sensitivity2(n,n);

    Eigen::Map<Eigen::VectorXd> map{grad.data(), long(grad.size())};
    if(doPositive) {
        sensitivity1.setFromTriplets(triplets1.begin(), triplets1.end());
        map -= l1.transpose()*sensitivity1;
    }
    if(doNegative) {
        sensitivity2.setFromTriplets(triplets2.begin(), triplets2.end());
        map -= l2.transpose()*sensitivity2;
    }

    Debug{} << "Gradient norm : " << map.norm();
    //Debug{} << "l norm : " << l.norm();
    //Debug{} << "Ap*x - bp norm : " << (Ap*x - bp).norm();
    //Debug{} << "Ap - Ap norm : " << (Ap - A).norm();
    //Debug{} << "bp - b norm : " << (bp - b).norm();
    //Debug{} << "grad x norm" << gradX.norm();
    //Debug{} << "x norm" << x.norm();
}

template<class Scalar>
void computeEnergy(
        DiffuseYamabe& yamabe,
        Eigen::Matrix<Scalar, Eigen::Dynamic, 1> const& x1,
        Eigen::Matrix<Scalar, Eigen::Dynamic, 1> const& x2,
        Eigen::Matrix<Scalar, Eigen::Dynamic, 1>& gradX1,
        Eigen::Matrix<Scalar, Eigen::Dynamic, 1>& gradX2,
        Scalar& out1,
        Scalar& out2) {

    using Vec3 = Math::Vector3<Scalar>;
    using VecX = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;

    Mesh& mesh = yamabe.mesh;

    bool doPositive = yamabe.positivePhase;
    bool doNegative = yamabe.negativePhase;
    bool computeDerivative = gradX1.size() > 0 || gradX2.size() > 0;

    if(yamabe.energy == EnergyType::Dirichlet) {
        for(Face face : mesh.faces()) {
            Math::Vector3<Scalar> grad1, grad2;
            Scalar rescaling = (yamabe.curvatureRescaling ? yamabe.getRescalingFactor(face) : 1.);

            for(HalfEdge he : face.halfEdges()) {
                Vertex v = he.next().tip();
                if(doPositive)
                    grad1 += x1[v.idx]*Vec3{mesh.gradient[he]};
                if(doNegative)
                    grad2 += x2[v.idx]*Vec3{mesh.gradient[he]};
            }

            if(computeDerivative) {
                for(HalfEdge he : face.halfEdges()) {
                    Vertex v = he.next().tip();
                    if(doPositive)
                        gradX1[v.idx] += rescaling*mesh.faceArea[face]*2*Math::dot(grad1, Vec3{mesh.gradient[he]});
                    if(doNegative)
                        gradX2[v.idx] += rescaling*mesh.faceArea[face]*2*Math::dot(grad2, Vec3{mesh.gradient[he]});
                }
            }

            Scalar grad1NormSquared = grad1.dot()*rescaling;
            Scalar grad2NormSquared = grad2.dot()*rescaling;

            out1 += face.area()*grad1NormSquared;
            out2 += face.area()*grad2NormSquared;
        }
    } else if(yamabe.energy == EnergyType::Hencky) {
        for(Vertex v : mesh.vertices()) {
            Scalar rescaling = (yamabe.curvatureRescaling ? yamabe.getRescalingFactor(v) : 1.);
            if(doPositive)
                out1 += x1[v.idx]*x1[v.idx]*mesh.integral[v];
            if(doNegative)
                out2 += x2[v.idx]*x2[v.idx]*mesh.integral[v];
            if(computeDerivative) {
                if(doPositive)
                    gradX1[v.idx] += 2*x1[v.idx]*mesh.integral[v]*rescaling;
                if(doNegative)
                    gradX2[v.idx] += 2*x2[v.idx]*mesh.integral[v]*rescaling;
            }
        }
    }
}

template<class Scalar>
void DiffuseYamabe::operator()(
         ArrayView<const Scalar> parameters, ArrayView<const Scalar> weights, Scalar& out,
         ArrayView<Scalar> gradP, ArrayView<Scalar> gradW) {

    using Vec3 = Math::Vector3<Scalar>;
    using VecX = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;

    size_t n = mesh.vertexCount();

    Eigen::SparseMatrix<Scalar> A1, A2;
    SolverType<Scalar> solver1, solver2;
    VecX b1, b2, x1, x2;
    Scalar lambda = lambdaWeight*(*scaling);
    assembleAndSolve<Scalar>(*this, lambda, parameters, weights, A1, solver1, A2, solver2, x1, x2, b1, b2);

    VecX gradX1{gradP ? n : 0};
    VecX gradX2{gradP ? n : 0};
    gradX1.setZero();
    gradX2.setZero();

    Scalar out1 = 0, out2 = 0;
    computeEnergy(*this, x1, x2, gradX1, gradX2, out1, out2);
    out += out1 + out2;

    if(gradP) {
        if constexpr (std::is_same_v<double, Scalar>) {
            shapeDerivative(*this, parameters, lambda, A1, solver1, A2, solver2, x1, x2, gradX1, gradX2, gradP);
        }
    }
}

size_t DiffuseYamabe::numParameters() const { return mesh.vertexCount(); }

void DiffuseYamabe::drawImGuiOptions(VisualizationProxy& proxy) {
    bool redraw = false;

    if(ImGui::Checkbox("Positive Phase", &positivePhase)) {
        redraw |= drawSolution && positivePhase;
    }

    ImGui::SameLine();

    if(ImGui::Checkbox("Negative Phase", &negativePhase)) {
        redraw |= drawSolution && negativePhase;
    }

    if(ImGui::Checkbox("Draw Solution", &drawSolution)) {
        if(drawSolution) {
            proxy.setCallbacks(
                    [this, &proxy](Node node) {
                        Eigen::SparseMatrix<double> A1, A2;
                        Eigen::VectorXd b1,b2,x1,x2;
                        SolverType<double> solver1, solver2;
                        double lambda = (*scaling)*lambdaWeight;
                        assembleAndSolve<double>(*this, lambda, node.phasefield(), node.temporary(), A1, solver1, A2, solver2, x1, x2, b1, b2);
                        Eigen::VectorXd x{mesh.vertexCount()};
                        x.setZero();
                        if(positivePhase) x += x1;
                        if(negativePhase) x += x2;
                        proxy.drawValuesNormalized({x.data(), size_t(x.size())});
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
                        constexpr double delta = 0.1;
                        auto phasefield = node.phasefield();

                        VertexData<size_t> map{DirectInit, phasefield.size(), Invalid};
                        VertexData<bool> inInterior{DirectInit, phasefield.size(), false};

                        size_t n = 0;
                        for(Vertex v : mesh.vertices()) {
                            if((positivePhase && phasefield[v] > 0) || (negativePhase && phasefield[v] < 0)) {
                                inInterior[v] = true;
                                map[v] = n++;
                            }
                        }

                        Array<Eigen::Triplet<double>> triplets;
                        Eigen::SparseMatrix<double> A(n,n);
                        Eigen::VectorXd b(n);

                        for(Vertex v : mesh.vertices()) {
                            size_t idx = map[v];

                            if(idx == Invalid) continue;
                            bool onBoundary = v.onBoundary();
                            if(!onBoundary) {
                                for(Vertex w : v.adjacentVertices()) {
                                    onBoundary |= map[w] == Invalid; /* if adjacent vertex is not in phase, the vertex is a boundary vertex */
                                }
                            }
                            inInterior[v] = !onBoundary;

                            if(!onBoundary) {
                                b[idx] = mesh.gaussianCurvature[v]*mesh.integral[v];
                            } else {
                                b[idx] = 0;
                                arrayAppend(triplets, InPlaceInit, int(idx), int(idx), 1.);
                            }
                        }


                        for(Face face : mesh.faces()) {
                            for(HalfEdge he1 : face.halfEdges()) {
                                Vertex v = he1.next().tip();
                                if(!inInterior[v]) continue;
                                for(HalfEdge he2 : face.halfEdges()) {
                                    Vertex w = he2.next().tip();
                                    if(!inInterior[w]) continue;
                                    double value = Math::dot(mesh.gradient[he2], mesh.gradient[he1])*face.area();
                                    //CORRADE_ASSERT(map[v] < n && map[w] < n, "Indices out of bound",);
                                    arrayAppend(triplets, InPlaceInit, map[v], map[w], value);
                                }
                            }
                        }

                        A.setFromTriplets(triplets.begin(), triplets.end());

#ifdef PHASEFIELD_WITH_SUITESPARSE
                        Eigen::UmfPackLU<Eigen::SparseMatrix<double>> solver;
#else
                        Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
#endif

                        solver.compute(A);
                        handleSolverInfo(solver.info());
                        Eigen::VectorXd x = solver.solve(b);

                        Debug{} << "Solution norm" << x.norm();

                        Array<double> solution{NoInit, phasefield.size()};
                        for(size_t i = 0; i < phasefield.size(); ++i) {
                            size_t idx = map[i];
                            if(idx != Invalid) {
                                solution[i] = x[idx];
                            } else {
                                solution[i] = 0;
                            }
                        }

                        proxy.drawValuesNormalized(solution);

                    },
                    [this]{ drawSolutionThresholded = false; });
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

    ImGui::SameLine();

    if(ImGui::Checkbox("Curvature Rescaling", &curvatureRescaling)) {
        redraw |= drawSolutionThresholded || drawSolution || drawSolutionGradient;
    }

    if(redraw) proxy.redraw();
}

double DiffuseYamabe::getRescalingFactor(Face f) const {
    double c = 0;
    constexpr double eps = 1e-4;
    for(Vertex v : f.vertices()) {
        c += mesh.gaussianCurvature[v];
    }
    c /= 3.;
    double r = 1./(Math::sqrt(Math::abs(c)) + eps);
    double a = r*r*r;
    assert(!std::isnan(a));
    return a;
}

double DiffuseYamabe::getRescalingFactor(Vertex v) const {
    double c = 0;
    constexpr double eps = 1e-4;
    double r = 1./(Math::sqrt(Math::abs(mesh.gaussianCurvature[v])) + eps);
    double a = r*r*r;
    assert(!std::isnan(a));
    return a;
}


DEFINE_FUNCTIONAL_CONSTRUCTOR(DiffuseYamabe)
DEFINE_FUNCTIONAL_OPERATOR(DiffuseYamabe, double)

#ifdef PHASEFIELD_WITH_ADOLC
DEFINE_FUNCTIONAL_OPERATOR(DiffuseYamabe, adouble)
#endif

}
