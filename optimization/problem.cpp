//
// Created by janos on 28.03.20.
//

#include "problem.hpp"
#include "normalizeInto.hpp"
#include "tag.h"

#include <Corrade/Utility/Assert.h>
#include <Corrade/Utility/Algorithms.h>
#include <Corrade/Containers/GrowableArray.h>

#include <adolc/adouble.h>
#include <adolc/taping.h>
#include <adolc/adolc_sparse.h>

using namespace Corrade;
using namespace Magnum;

namespace solver {

namespace {
void add(SparseMatrix const& a, SparseMatrix& b){
    std::size_t idx = 0;
    for(std::size_t i = 0; i < a.nnz; ++i){
        UnsignedInt row = a.rows[i];
        UnsignedInt col = a.cols[i];
        while(row != b.rows[idx] || col != b.cols[idx])
            ++idx;
        b.values[idx] += a.values[i];
    }
}

void append(SparseMatrix& a, SparseMatrix const& b, UnsignedInt rowOffset){

    Containers::arrayResize(a.values, a.nnz + b.nnz);
    Containers::arrayResize(a.rows, a.nnz + b.nnz);
    Containers::arrayResize(a.cols, a.nnz + b.nnz);

    for(std::size_t i = 0; i < b.nnz; ++i){
        a.values[i + a.nnz] = b.values[i];
        a.rows[i + a.nnz] = rowOffset + b.rows[i];
        a.cols[i + a.nnz] = b.cols[i];
    }

    a.nnz += b.nnz;
}
}


Problem::Problem() : tagL(getTag()), tagJ(getTag()), tagG(getTag()) {}
Problem::~Problem(){
    deleteTag(tagL);
    deleteTag(tagJ);
    deleteTag(tagG);
}

void Problem::evaluate(double const* parameters,
                       double* residual,
                       double* constr,
                       double* g,
                       SparseMatrix* j,
                       SparseMatrix* h,
                       double objectiveScale,
                       double const* lambdas
                       ) const {

    CORRADE_ASSERT(!objectives.empty(), "Problem : need at least on objective function in problem", );

    auto n = numParameters();
    auto m = numParameters();

    if(g){
        std::fill_n(g, n, 0.);
    }

    Containers::Array<double> grad(Containers::ValueInit, numParameters());
    for(std::size_t i = 0; i < objectives.size(); ++i){
        double r;
        double* gradient = g ? grad.data() : nullptr;
        objectives[i](parameters, &r, gradient, nullptr);
        if(g){
            for(int k = 0; k < n; ++k){
                g[k] += gradient[i];
            }
        }
    }

    if(j){
        UnsignedInt* rInd = const_cast<UnsignedInt*>(hessian.rows.data());
        UnsignedInt* cInd = const_cast<UnsignedInt*>(hessian.cols.data());
        int nnz = hessian.nnz;
        double* val = const_cast<double*>(hessian.values.data());
        int options[4] = {0,0,0,0}; /* all default */
        sparse_jac(tagJ, m, n, 1, const_cast<double*>(parameters), &nnz, &rInd, &cInd, &val, options);
    }

    if(h){
        Containers::Array<double> parametricValues(m + 1);
        parametricValues[0] = objectiveScale;
        for(std::size_t i = 1; i < m + 1; ++i){
            parametricValues[i] = lambdas[i];
        }
        set_param_vec(tagL, m + 1, parametricValues);

        UnsignedInt* rInd = const_cast<UnsignedInt*>(hessian.rows.data());
        UnsignedInt* cInd = const_cast<UnsignedInt*>(hessian.cols.data());
        int nnz = hessian.nnz;
        double* val = const_cast<double*>(hessian.values.data());
        int options[2] = {0,0}; /* all default */
        sparse_hess(tagL, n, 1, const_cast<double*>(parameters), &nnz, &rInd, &cInd, &val, options);
    }

}

std::size_t Problem::numParameters() const {
    if (objectives.empty() && constraints.empty())
        return 0;
    auto numParams = objectives.front().numParameters();
#ifdef CORRADE_ASSERT
    for (auto const& f : objectives)
        CORRADE_INTERNAL_ASSERT(f.numParameters() == numParams);
    for (auto const& c : constraints)
        CORRADE_INTERNAL_ASSERT(c.numParameters() == numParams);
#endif
    return numParams;
}

std::size_t Problem::numConstraints() const {
    std::size_t m = 0;
    for(auto const& c : constraints)
        m += c.numResiduals();
    return m;
}

void Problem::updateInternalDataStructures(Containers::Array<double> const& x) {

    for(auto& f : objectives)
        f.updateInternalDataStructures();
    for(auto& f : objectives)
        f.updateInternalDataStructures();

    auto n = numParameters();
    auto m = numConstraints();
    Containers::Array<adouble> cs(m);
    Containers::Array<adouble> xa(n);
    double dummy;

    UnsignedInt *rInd, *cInd;
    double* val;

    {
        trace_on(tagG);

        for(std::size_t idx = 0; idx < n; idx++)
            xa[idx] <<= x[idx];

        adouble objValue = 0;
        for(Functional const& f : objectives){
            adouble fr;
            f(xa, &fr);
            objValue += fr;
        }

        objValue >>= dummy;

        trace_off();
    }

    {
        trace_on(tagJ);

        for(std::size_t idx = 0; idx < n; idx++)
            xa[idx] <<= x[idx];

        int offset = 0;
        for(auto const& f : constraints){
            auto fm = f.numResiduals();
            f(xa, cs.begin() + offset);
            offset += fm;
        }

        for(std::size_t i = 0; i < m; ++i)
            cs[i] >>= dummy;

        trace_off();

        int options[4] = {0,0,0,0}; /* all default */
        sparse_jac(tagJ, m, n, 0, x.data(), &jacobian.nnz, &rInd, &cInd, &val, options);
        jacobian.values = Containers::Array(val, jacobian.nnz);
        jacobian.rows = Containers::Array(rInd, jacobian.nnz);
        jacobian.cols = Containers::Array(cInd, jacobian.nnz);
    }

    {
        trace_on(tagL);

        for(std::size_t idx = 0; idx < n; idx++)
            xa[idx] <<= x[idx];

        adouble objValue = 0;

        for(Functional const& f : objectives){
            adouble fr;
            f(xa, &fr);
            objValue += fr;
        }
        objValue *= mkparam(1.);

        int offset = 0;
        for(std::size_t i = 0; i < constraints.size(); ++i){
            auto fm = constraints[i].numResiduals();
            objectives[i](xa, cs.begin() + offset);
            for(std::size_t j = 0; j < fm; ++j){
                objValue += cs[offset + j] * mkparam(1.);
            }
            offset += fm;
        }
        objValue >>= dummy;

        trace_off();

        int options[2] = {0,0}; /* all default */
        sparse_hess(tagL, n, 0, x.data(), &hessian.nnz, &rInd, &cInd, &val, options);
        hessian.values = Containers::Array(val, hessian.nnz);
        hessian.rows = Containers::Array(rInd, hessian.nnz);
        hessian.cols = Containers::Array(cInd, hessian.nnz);
    }
}

}
