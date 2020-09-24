//
// Created by janos on 20.04.20.

#include "LossFunctions.h"
#include "SmartEnum.h"

#include <Corrade/Utility/StlMath.h>
#include <Corrade/Utility/Assert.h>

//#include <adolc/adouble.h>

#include <imgui.h>
#include <utility>
#include <new>

namespace Phasefield {

template<class T>
LossFunction::LossFunction(T f) {
    erased = ::operator new(sizeof(T), std::align_val_t(alignof(T)));
    ::new(erased) T(std::move(f));
    //ad = +[](void* e, adouble const& x, adouble& y) {
    //    (*static_cast<T*>(e))(x, y);
    //};
    loss = +[](void* e, double x, double ys[3]) {
        (*static_cast<T*>(e))(x, ys);
    };
    destroy = +[](void* e) {
        static_cast<T*>(e)->~T();
        ::operator delete(e, sizeof(T), std::align_val_t(alignof(T)));
    };

    if constexpr ( requires {T::type(); })
        lossType = T::type();
    else
        lossType = LossFunctionType::Unknown;

    /* optional */
    //constexpr bool hasSettings = is_valid(f)([](auto&& p) constexpr -> decltype(f.drawSettings()){});
    if constexpr(requires { f.drawSettings(); }) {
        draw = +[](void* e) { static_cast<T*>(e)->drawSettings(); };
    }
}

LossFunction::LossFunction(LossFunctionType::Value t) {
    switch(t) {
        case LossFunctionType::Cauchy : *this = CauchyLoss{}; break;
        case LossFunctionType::Trivial : *this = TrivialLoss{}; break;
        case LossFunctionType::Quadratic : *this = QuadraticLoss{}; break;
        default : CORRADE_ASSERT(false, "Unknown loss type",);
    }
}

void LossFunction::drawSettings() {
    if(ImGui::BeginCombo("##loss", LossFunctionType::to_string(lossType))) {
        for(auto l : LossFunctionType::range) {
            if(l == LossFunctionType::Unknown) continue;
            bool isSelected = l == lossType;
            if(ImGui::Selectable(LossFunctionType::to_string(l), isSelected))
                *this = LossFunction(l);
            if(isSelected) {
                ImGui::SetItemDefaultFocus();
            }
        }

        ImGui::EndCombo();
    }

    if(draw) {
        draw(erased);
    }
    ImGui::InputDouble("Weight", &weight);
}

LossFunction::~LossFunction() {
    if(erased) {
        destroy(erased);
    }
}

LossFunction::LossFunction(LossFunction&& other) noexcept {
    swap(other);
}

LossFunction& LossFunction::operator=(LossFunction&& other) noexcept {
    swap(other);
    return *this;
}

void LossFunction::swap(LossFunction& other) {
    //std::swap(ad, other.ad);
    std::swap(loss, other.loss);
    std::swap(destroy, other.destroy);
    std::swap(erased, other.erased);
    std::swap(lossType, other.lossType);
    std::swap(weight, other.weight);
    std::swap(draw, other.draw);
}

void LossFunction::operator()(double const& in, double out[3]) const {
    loss(erased, in, out);
    for(size_t i = 0; i < 3; ++i) out[i] *= weight;
}

//void LossFunction::operator()(adouble const& x, adouble& y) const {
//    ad(erased, x, y);
//    y *= weight;
//}

void TrivialLoss::operator()(double const& in, double out[3]) const {
    out[0] = in;
    out[1] = 1.;
    out[2] = 0.;
}

void QuadraticLoss::operator()(double const& in, double out[3]) const {
    out[0] = in*in;
    out[1] = 2*in;
    out[2] = 2.;
}

void CauchyLoss::operator()(double const& in, double out[3]) const {
    double sum = in + 1.;
    double inv = 1./sum;
    out[0] = log(sum);
    out[1] = inv;
    out[2] = -1.*inv*inv;
}


/* explicit instantiations */
template LossFunction::LossFunction(TrivialLoss);

template LossFunction::LossFunction(CauchyLoss);

template LossFunction::LossFunction(QuadraticLoss);

}