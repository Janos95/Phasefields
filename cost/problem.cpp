//
// Created by janos on 28.03.20.
//

#include "problem.hpp"
#include "normalizeInto.hpp"

#include <Corrade/Utility/Assert.h>
#include <Corrade/Utility/Algorithms.h>
#include <Corrade/Containers/Array.h>

#include <mutex>
#include <tbb/task_group.h>
#include <Magnum/Trade/MeshData.h>

using namespace Corrade;
using namespace Magnum;

namespace solver {


bool Problem::evaluate(double const *parameters,
                       double *costF,
                       double *jacF,
                       double *costC,
                       double *jacC) const {
    CORRADE_ASSERT(!functionals.empty(), "Problem : no functionals in problem", false);
    tbb::task_group g;

    if(costC) *costC = 0.;
    if(costF) *costF = 0.;

    if (jacF)
        std::fill_n(jacF, numParameters(), 0.);

    if(jacC)
        std::fill_n(jacC, numParameters() * constraints.size(), 0.);

    std::mutex mCostF;
    std::mutex mJacF;
    std::mutex mCostC;
    std::mutex mJacC;

    auto eval = [](
            FunctionalD const& f,
            double const * const parameters,
            double* cost, std::mutex& mCost,
            double* jacobian, std::mutex& mJac)
    {
        if(!cost && !jacobian) return;
        auto n = f.numParameters();
        Containers::Array<Double> grad(Containers::NoInit, jacobian ? n : 0);
        Double residual = 0;
        f.evaluate(parameters, &residual, jacobian ? grad.data() : nullptr);

        Double rho[3], scaling[3] = {residual, 1., 0.};
        if(f.metaData->scaling){
            auto s = *f.metaData->scaling;
            for (auto& r : scaling) r *= s;
        }
        f.metaData->loss->Evaluate(scaling[0], rho);
        if (jacobian) {
            auto derivative = rho[1] * scaling[1];
            for(auto& x : grad) x *= derivative;
            f.metaData->visualizeGradient(grad);
            std::lock_guard l(mJac);
            for (int i = 0; i < n; ++i)
                jacobian[i] += grad[i];
        }
        if(cost){
            std::lock_guard l(mCost);
            *cost += rho[0];
        }
    };

    for (auto const& f: functionals)
        g.run([&]{ eval(*f, parameters, costF, mCostF, jacF, mJacF); });

    for (std::size_t i = 0; i < constraints.size(); ++i) {
        auto jac = jacC ? jacC + i * numParameters() : jacC; // dense jacobian
        g.run([&]{ eval(*constraints[i], parameters, costC, mCostC, jac, mJacC); });
    }

    g.wait();
    {
        std::size_t n = numParameters();
        if(!meshData || !jacF || !mutex) return true;
        std::lock_guard l(*mutex);
        auto coords = meshData->mutableAttribute(Mg::Trade::MeshAttribute::TextureCoordinates);
        auto xcoords = Cr::Containers::arrayCast<2, Mg::Float>(coords).slice<1>();
        if(flags & VisualizationFlag::Gradient){
            normalizeInto({jacF, n}, xcoords);
            *update |= VisualizationFlag::Gradient;
        }
        if(flags & VisualizationFlag::Phasefield){
            normalizeInto({parameters, n}, xcoords, -1., 1.);
            *update |= VisualizationFlag::Phasefield;
        }
    }
    return true;
}

int Problem::numParameters() const {
    if (functionals.empty() && constraints.empty())
        return 0;
    auto numParams = functionals.front()->numParameters();
#ifndef NDEBUG
    for (auto const& f : functionals)
        CORRADE_INTERNAL_ASSERT(f->numParameters() == numParams);
    for (auto const& c : constraints)
        CORRADE_INTERNAL_ASSERT(c->numParameters() == numParams);
#endif
    return numParams;
}

int Problem::numConstraints() const {
    return constraints.size();
}

void Problem::numericalGradient(double* parameters, double* jacF, double* jacC) {
    std::size_t n = numParameters();
    std::size_t m = numConstraints();
    Cr::Containers::Array<double> gradF(Cr::Containers::NoInit, n);
    Cr::Containers::Array<double> gradC(Cr::Containers::NoInit, m*n);
    Cr::Containers::Array<double> cs1(Cr::Containers::NoInit, m), cs2(Cr::Containers::NoInit, m);
    for (std::size_t i = 0; i < n; ++i) {
        double r1, r2;
        double p = parameters[i];
        parameters[i] += h;
        evaluate(parameters, &r2, nullptr, cs1.data(), nullptr);
        parameters[i] -= 2.*h;
        evaluate(parameters, &r1, nullptr, cs2.data(), nullptr);
        jacF[i] = (r2 - r1) / (2.*h);
        for (std::size_t j = 0; j < m; ++j) {
            jacC[i + n * j] = (cs1[j] - cs2[j]) / (2.*h);
        }
        //restore old value
        parameters[i] = p;
    }
}

}
