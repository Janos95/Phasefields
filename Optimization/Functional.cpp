//
// Created by janos on 6/25/20.
//


#include "Functional.h"
#include "../Cost/SparseMatrix.h"
#include "Tag.h"

#include <Corrade/Containers/GrowableArray.h>
#include <Magnum/Math/Functions.h>

#include <adolc/adouble.h>
#include <adolc/taping.h>
#include <adolc/adolc_sparse.h>

//#include <imgui.h>

namespace Phasefield {

using namespace Magnum;

Functional::~Functional() {
    deleteTag(tag);
    if(erased) {
        destroy(erased);
        //deleteTag(tag);
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

void Functional::operator()(Containers::ArrayView<const double> parameters,
                            Containers::ArrayView<const double> weights,
                            double& out,
                            Containers::ArrayView<double> gradP,
                            Containers::ArrayView<double> gradW) const {

    auto n = numParameters();
    double cost = 0;
    evalWithGrad(erased, parameters, weights, cost, gradP, gradW);

    /* apply scaling and loss function */
    Mg::Double rho[3], phi[3] = {cost, 1., 0.};
    if(scaling) {
        auto s = *scaling;
        for(auto& r : phi) r *= s;
    }

    loss(phi[0], rho);

    out += rho[0];

    if(gradP || gradW) {
        double lossGrad = rho[1]*phi[1];
        for(std::size_t i = 0; i < n; ++i) {
            if(gradP) gradP[i] *= lossGrad;
            if(gradW) gradW[i] *= lossGrad;
        }
    }
}

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

bool Functional::drawImGuiOptions(VisualizationProxy& proxy) {
    if(options) {
        //ImGui::Checkbox("Check Derivatives", &checkDerivatives);
        return options(erased, proxy);
    }
    return false;

}

void Functional::updateInternalDataStructures() {
    update(erased);
    //isFirstEvaluation = true;
}

[[nodiscard]] std::size_t Functional::numParameters() const {
    return params(erased);
}

void swap(Functional& f1, Functional& f2) {
    std::swap(f1.erased, f2.erased);
    std::swap(f1.destroy, f2.destroy);
    std::swap(f1.params, f2.params);
    //std::swap(f1.residuals, f2.residuals);
    std::swap(f1.vis, f2.vis);
    std::swap(f1.off, f2.off);
    std::swap(f1.update, f2.update);
    std::swap(f1.options, f2.options);
    std::swap(f1.evalWithGrad, f2.evalWithGrad);
    //std::swap(f1.ad, f2.ad);
    //std::swap(f1.tag, f2.tag);
    //std::swap(f1.isFirstEvaluation, f2.isFirstEvaluation);
    std::swap(f1.alwaysRetape, f2.alwaysRetape);

    /* two phase lookup */
    using std::swap;
    swap(f1.loss, f2.loss);
    swap(f1.scaling, f2.scaling);
}

Functional::Functional(Functional&& other) noexcept {
    using std::swap;
    swap(*this, other);
}

Functional& Functional::operator=(Functional&& other) noexcept {
    using std::swap;
    swap(*this, other);
    return *this;
}

/* this is called from the gui thread so we can update some opengl stuff if we want to */
void Functional::updateVisualization(VisualizationProxy& proxy) {
    if(vis)
        vis(erased, proxy);
}

void Functional::turnVisualizationOff() {
    if(off)
        off(erased);
}

}
