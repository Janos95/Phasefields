//
// Created by janos on 28.03.20.
//

#include "problem.hpp"
#include "normalizeInto.hpp"

#include <Corrade/Utility/Assert.h>
#include <Corrade/Containers/Array.h>

#include <mutex>
#include <tbb/task_group.h>
#include <Magnum/Trade/MeshData.h>

using namespace Corrade;
using namespace Magnum;

namespace solver {

bool Problem::evaluate(double const *parameters,
                       double *cost,
                       double *jacobians) const {
    CORRADE_ASSERT(!functionals.empty(), "Problem : no functionals in problem", false);
    const std::size_t n = numParameters();

    *cost = 0.;
    if (jacobians)
        std::fill_n(jacobians, numParameters(), 0.);

    std::mutex m1;
    std::mutex m2;

    //tbb::task_group g;
    for (auto const& functional: functionals) {
        //g.run([&] {
            Containers::Array<Double> grad(Containers::NoInit, jacobians ? n : 0);
            Double residual = 0;
            functional->evaluate(parameters, &residual, jacobians ? grad.data() : nullptr);

            Double rho[3], scaling[3] = {residual, 1., 0.};
            if(functional->metaData->scaling){
                auto s = *functional->metaData->scaling;
                for (auto& r : scaling) r *= s;
            }
            functional->metaData->loss->Evaluate(scaling[0], rho);
            if (jacobians) {
                auto derivative = rho[1] * scaling[1];
                for(auto& x : grad) x *= derivative;
                functional->metaData->visualizeGradient(grad);
                std::lock_guard l(m1);
                for (int i = 0; i < n; ++i)
                    jacobians[i] += grad[i];
            }
            std::lock_guard l(m2);
            *cost += rho[0];
       // });
    }
    //g.wait();
    {
        if(!meshData || !jacobians || !mutex) return true;
        std::lock_guard l(*mutex);
        auto coords = meshData->mutableAttribute(Mg::Trade::MeshAttribute::TextureCoordinates);
        auto xcoords = Cr::Containers::arrayCast<2, Mg::Float>(coords).slice<1>();
        if(flags & VisualizationFlag::Gradient){
            normalizeInto({jacobians, n}, xcoords);
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

}
