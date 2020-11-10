//
// Created by janos on 6/25/20.
//


#include "Functional.h"
#include "SparseMatrix.h"
#include "Tag.h"
#include "VisualizationProxy.h"
#include "Viewer.h"
#include "ScopedTimer/ScopedTimer.h"

#include <Corrade/Utility/ConfigurationGroup.h>

#include <Corrade/Containers/GrowableArray.h>
#include <Corrade/Utility/Algorithms.h>
#include <Corrade/Utility/FormatStl.h>
#include <Magnum/Math/Functions.h>
#include <Magnum/Math/FunctionsBatch.h>

#ifdef PHASEFIELD_WITH_ADOLC
#include <adolc/adouble.h>
#include <adolc/taping.h>
#include <adolc/drivers/drivers.h>
#endif

#include <imgui.h>

namespace Phasefield {

Functional::~Functional() {
    //deleteTag(tag);
    if(erased) {
        destroy(erased);
        deleteTag(tag);
    }
}

void Functional::operator()(ArrayView<const double> parameters,
                            ArrayView<const double> weights,
                            double& out,
                            ArrayView<double> gradP,
                            ArrayView<double> gradW) const {

    ScopedTimer timer{FunctionalType::to_string(functionalType), false};


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

#ifdef PHASEFIELD_WITH_ADOLC
    if(checkDerivatives && gradP) {
        Debug{} << "Checking Derivatives";

        CORRADE_INTERNAL_ASSERT(tag != Invalid);
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
    }
#endif

#ifndef NDEBUG
    if(!check({&out, 1})) {
        Debug{} << "NaN in cost of functional" << FunctionalType::to_string(functionalType);
    }
    if(gradP) {
        if(!check(gradP)) {
            Debug{} << "NaN in derivative of functional" << FunctionalType::to_string(functionalType);
        }
    }
#endif

}

void Functional::drawImGuiOptions(VisualizationProxy& proxy) {

    ImGui::Text(FunctionalType::to_string(functionalType));

    //DrawableType type;
    loss.drawSettings();

    if(ImGui::Checkbox("Draw Derivative", &drawGradient)) {
        if(drawGradient) {
            if(functionalType == FunctionalType::AreaRegularizer) {
                proxy.viewer.totalArea = proxy.viewer.currentNode.integrateWeight(proxy.viewer.mesh);
            }
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
    std::swap(scaling, other.scaling);
    std::swap(load, other.load);
    std::swap(save, other.save);

    loss.swap(other.loss);
}

Functional::Functional(Functional&& other) noexcept {
    swap(other);
}

Functional& Functional::operator=(Functional&& other) noexcept {
    swap(other);
    return *this;
}

bool Functional::check(ArrayView<const double> data) const {
    for(double x : data) {
        if(std::isnan(x)) return false;
    }
    return true;
}

void Functional::loadParameters(Cr::Utility::ConfigurationGroup const& group) {
    if(group.hasValue("weight")) {
        loss.weight = group.value<double>("weight");
    }
    if(load) {
        load(erased, group);
    }
}

void Functional::saveParameters(Cr::Utility::ConfigurationGroup& group) const {
    if(save) {
        group.setValue("weight", loss.weight);
        save(erased, group);
    }
}

}
