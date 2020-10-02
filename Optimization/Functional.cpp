//
// Created by janos on 6/25/20.
//


#include "Functional.h"
#include "SparseMatrix.h"
#include "Tag.h"
#include "VisualizationProxy.h"
#include "Viewer.h"

#include <Corrade/Containers/GrowableArray.h>
#include <Corrade/Utility/Algorithms.h>
#include <Corrade/Utility/FormatStl.h>
#include <Magnum/Math/Functions.h>
#include <Magnum/Math/FunctionsBatch.h>

#include <adolc/adouble.h>
#include <adolc/taping.h>
#include <adolc/drivers/drivers.h>
//#include <adolc/adolc_sparse.h>



#include <imgui.h>

namespace Phasefield {

Functional::~Functional() {
    //deleteTag(tag);
    if(erased) {
        destroy(erased);
        deleteTag(tag);
    }
}

//void Functional::operator()(const double* parameters, const double* weights, double* out, double* gradP,
//                            double* gradW) const {
//    eval(erased, parameters, weights, out, gradP, gradW);
//}
//
//void Functional::operator()(const adouble* params, const adouble* weights, adouble* out) const {
//    ad(erased, params, weights, out);
//}


#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60



void printProgress(double percentage) {
    int val = (int) (percentage * 100);
    int lpad = (int) (percentage * PBWIDTH);
    int rpad = PBWIDTH - lpad;
    printf("\r%3d%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
    fflush(stdout);
}

#pragma clang diagnostic push
#pragma ide diagnostic ignored "openmp-use-default-none"
void Functional::operator()(ArrayView<const double> parameters,
                            ArrayView<const double> weights,
                            double& out,
                            ArrayView<double> gradP,
                            ArrayView<double> gradW) const {

    size_t n = numParameters();

    Array<double> gradPWithoutLoss{gradP ? n : 0};
    Array<double> gradWWithoutLoss{gradW ? n : 0};

    double cost = 0;
    evalWithGrad(erased, parameters, weights, cost, gradPWithoutLoss, gradWWithoutLoss);

    /* apply scaling and loss function */
    double rho[3], phi[3] = {cost, 1., 0.};
    if(scaling) {
        double s = *scaling;
        for(double& r : phi) r *= s;
    }

    loss(phi[0], rho);

    out += rho[0];

    double lossGrad = rho[1]*phi[1];
    if(gradP || gradW) {
        for(size_t i = 0; i < n; ++i) {
            if(gradP) gradP[i] += lossGrad*gradPWithoutLoss[i];
            if(gradW) gradW[i] += lossGrad*gradWWithoutLoss[i];
        }
    }

    if(checkDerivatives && gradP) {
        Debug{} << "Checking Derivatives";

        CORRADE_INTERNAL_ASSERT(tag != -1);
        trace_on(tag);

        Array<adouble> paramsAd{n};
        Array<adouble> weightsAd{n};
        adouble residual = 0;

        Array<double> gradientsAd(2*n);
        Array<double> variables(2*n);

        for(size_t i = 0; i < n; i++) {
            paramsAd[i] <<= parameters[i];
            variables[i] = parameters[i];
        }
        for(size_t i = 0; i < n; i++) {
            weightsAd[i] <<= weights[i];
            variables[i + n] = weights[i];
        }
        ad(erased, paramsAd, weightsAd, residual);
        if(scaling) residual *= *scaling;
        loss(residual, residual);

        double dummy;
        residual >>= dummy;

        trace_off();

        gradient(tag,2*n,variables.data(),gradientsAd.data());

        auto gradParamsAd = gradientsAd.prefix(n);

        double tol = 1e-6;
        for(size_t i = 0; i < n; ++i) {
            if(Math::abs(gradParamsAd[i] - lossGrad*gradPWithoutLoss[i]) > tol) {
                Debug{} << "Gradient Error in" << FunctionalType::to_string(functionalType);
                Debug{} << Cr::Utility::formatString("(Analytic = {}, AD = {}) at index {}\n", gradPWithoutLoss[i], gradParamsAd[i], i).c_str();
                CORRADE_INTERNAL_ASSERT(false);
            }
        }

        //Array<double> perturbedParams{NoInit, n};
        //Cr::Utility::copy(parameters, perturbedParams);

        //double h = 1e-6;
        //double tol = 1e-4;

        //std::atomic_size_t counter = 0;

//#pragm//a omp parallel for
        //for(size_t i = 0; i < n; ++i) {
        //    double f1 = 0, f2 = 0;
        //    double old = perturbedParams[i];
        //    perturbedParams[i] += h;
        //    evalWithGrad(erased, perturbedParams, weights, f2, nullptr, nullptr);
        //    perturbedParams[i] = old - h;
        //    evalWithGrad(erased, perturbedParams, weights, f1, nullptr, nullptr);
        //    perturbedParams[i] = old;

        //    double gradNumeric = (f2 - f1)/(2*h);
        //    double gradAnalytic = gradPWithoutLoss[i];

        //    if(Math::abs(gradNumeric - gradPWithoutLoss[i]) > tol) {
        //        Debug{} << "Gradient Error in" << FunctionalType::to_string(functionalType);
        //        Debug{} << Cr::Utility::formatString("(Analytic = {}, Numeric = {}) at index {}\n", gradAnalytic, gradNumeric, i).c_str();
        //        CORRADE_INTERNAL_ASSERT(false);
        //    }

        //    double progress = double(counter++)/double(n);
        //}
    }
}
#pragma clang diagnostic pop

//void Functional::operator()(
//        double const* parameter,
//        double* residuals,
//        double* gradient,
//        SparseMatrix* hessian) const {
//
//    //CORRADE_ASSERT(!hessian, "Functional : hessian not yet supported", );
//    //CORRADE_ASSERT(numResiduals() == 1, "Functional : can only use this overload if you have exactly one output residual",);
//
//    //bool needsTapingForHessian = alwaysRetape && hessian && nullptr == evalWithHessian;
//    //bool needsTapingForJacobian = alwaysRetape && jacobian && nullptr == evalWithJacobian;
//    //bool needsTapingForCheck = alwaysRetape && checkDerivatives;
//    //bool needsTaping = needsTapingForHessian || needsTapingForJacobian;
//    //int repeat = 1;
//
//    //if(needsTaping || needsTapingForCheck || isFirstEvaluation){
//    //    isFirstEvaluation = false;
//    //    tapeEvaluation(parameter);
//    //    repeat = 0;
//    //}
//
//    //std::size_t m = numParameters();
//    //std::size_t n = numResiduals();
//
//    //int jacobianOptions[4] = {0,0,0,0};
//    //int hessianOptions[2] = {0,0};
//
//    //double *hessianValues, *jacobianValues;
//    //unsigned int *hessianRowIndices, *hessianColumnIndices;
//    //unsigned int *jacobianRowIndices, *jacobianColumnIndices;
//    //int nnz;
//
//    //if(evalWithHessian) {
//    //    evalWithHessian(erased, parameter, residuals, jacobian, hessian);
//    //} else if(evalWithJacobian){
//    //    evalWithJacobian(erased, parameter, residuals, jacobian);
//    //    if(hessian) {
//    //        nnz = hessian->nnz;
//    //        set_param_vec(tag, m, const_cast<double*>(lambdas));
//    //        sparse_hess(tag, n, repeat, parameter, &nnz, &hessianRowIndices, &hessianColumnIndices, &hessianValues, hessianOptions);
//    //        if(!repeat){
//    //            hessian->nnz = nnz;
//    //            hessian->cols = Containers::Array(hessianColumnIndices, nnz);
//    //            hessian->rows = Containers::Array(hessianRowIndices, nnz);
//    //            hessian->values = Containers::Array(hessianValues, nnz);
//    //        }
//    //    }
//    //}
//
//    //if(evalWithJacobian){
//    //    evalWithJacobian(erased, parameter, residuals, jacobian);
//    //    if(hessian){
//    //        sparse_hess(tag, m, n, repeat, parameter, &nnz_jac, &rind_g, &cind_g, &jacval, options_g);
//    //    }
//    //    if(checkDerivatives && jacobian){
//    //        sparse_grad(tag, m, n, repeat, parameter, &nnz_jac, &rind_g, &cind_g, &jacval, options_g);
//    //        CORRADE_ASSERT();
//    //    }
//    //} else {
//    //    if(jacobian){
//    //        sparse_grad(tag, m, n, 1, x, &nnz_jac, &rind_g, &cind_g, &jacval, options_g);
//    //    }
//
//    //}
//
//    ///* apply scaling and loss function */
//    //for (std::size_t i = 0; i < m; ++i) {
//    //    Mg::Double rho[3], phi[3] = {residuals[i], 1., 0.};
//    //    if(scaling){
//    //        auto s = *scaling;
//    //        for (auto& r : phi) r *= s;
//    //    }
//
//    //    loss(phi[0], rho);
//
//    //    residuals[i] = rho[0];
//
//    //    if (gradient) {
//    //        double lossGrad = rho[1] * phi[1];
//    //        for(auto& x : jacobian->row(i))
//    //            x *= lossGrad;
//    //    }
//    //}
//}
//
//void Functional::tapeEvaluation(double const* x) const {
//    double dummy;
//    std::size_t m = numParameters();
//    std::size_t n = numResiduals();
//    trace_on(tag);
//    Containers::Array<adouble> xs(m);
//    Containers::Array<adouble> ys(m);
//    for (UnsignedInt i = 0; i < n; ++i)
//        xs[i] <<= x[i];
//    ad(erased, xs.data(), ys.data());
//    adouble y = 0.;
//    for(uint32_t i = 0; i < ys.size(); ++i){
//        ys[i] *= *scaling;
//        loss(ys[i], ys[i]);
//    }
//    y >>= dummy;
//    trace_off();
//}

void Functional::drawImGuiOptions(VisualizationProxy& proxy) {

    ImGui::Text(FunctionalType::to_string(functionalType));

    //DrawableType type;
    loss.drawSettings();

    if(ImGui::Checkbox("Draw Derivative", &drawGradient)) {
        if(drawGradient) {
            proxy.setCallbacks([this, &proxy](Node node) {
                auto parameters = node.phasefield();
                auto weights = node.temporary();
                Array<double> gradP(parameters.size());
                double cost = 0;
                (*this)(parameters, weights, cost, gradP, nullptr);
                proxy.drawValuesNormalized(gradP);
            },
            [this]{ drawGradient = false; });
        } else {
            proxy.setDefaultCallback();
        }
        proxy.redraw();
    }

    ImGui::SameLine();
    ImGui::Checkbox("Check Derivatives Using AD", &checkDerivatives);

    if(options) {
        options(erased, proxy);
    }

}

//void Functional::updateInternalDataStructures() {
//    update(erased);
//    //isFirstEvaluation = true;
//}

[[nodiscard]] std::size_t Functional::numParameters() const {
    return params(erased);
}

void Functional::swap(Functional& other) {
    std::swap(erased, other.erased);
    std::swap(destroy, other.destroy);
    std::swap(params, other.params);
    std::swap(options, other.options);
    std::swap(evalWithGrad, other.evalWithGrad);
    std::swap(functionalType, other.functionalType);
    std::swap(tag, other.tag);
    std::swap(ad, other.ad);

    loss.swap(other.loss);
    scaling.swap(other.scaling);
}

Functional::Functional(Functional&& other) noexcept {
    swap(other);
}

Functional& Functional::operator=(Functional&& other) noexcept {
    swap(other);
    return *this;
}

/* this is called from the gui thread so we can update some opengl stuff if we want to */
//void Functional::updateVisualization(VisualizationProxy& proxy) {
//    if(vis)
//        vis(erased, proxy);
//}

//void Functional::turnVisualizationOff() {
//    if(off)
//        off(erased);
//}

}
