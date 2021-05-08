//
// Created by janos on 11/3/20.
//

#pragma once

#include "CudaCG.h"
#include "Types.h"

#include <Corrade/Utility/Debug.h>

#include <Eigen/SparseCore>
#include <Eigen/SparseLU>
#include <Eigen/IterativeLinearSolvers>
//#include <Eigen/PardisoSupport>

#ifdef PHASEFIELD_WITH_SUITESPARSE
#include <Eigen/UmfPackSupport>
#include <Eigen/CholmodSupport>
#endif

#ifdef PHASEFIELD_WITH_ADOLC
#include <adolc/adouble.h>

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

#endif

namespace Phasefield {

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

template<class Scalar>
struct SelectSolver {
#ifndef MAGNUM_TARGET_WEBGL
#ifdef PHASEFIELD_WITH_SUITESPARSE
    //using type = Eigen::UmfPackLU<Eigen::SparseMatrix<Scalar>>;
    //using type = Eigen::CholmodSimplicialLDLT<Eigen::SparseMatrix<Scalar>>;
    //using type = Eigen::PardisoLDLT<Eigen::SparseMatrix<Scalar>, Eigen::Upper>;
    //using type = CudaCG;
    //using type = Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper>;
    using type = Eigen::CholmodSupernodalLLT<Eigen::SparseMatrix<Scalar>>;
#else
    using type = Eigen::SparseLU<Eigen::SparseMatrix<Scalar>>;
#endif
#else
    using type = Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper>;
#endif
};

#ifdef PHASEFIELD_WITH_ADOLC
template<>
struct SelectSolver<adouble> {
    using type = Eigen::SparseLU<Eigen::SparseMatrix<adouble>>;
};
#endif

template<class Scalar>
using SolverType = typename SelectSolver<Scalar>::type;

}