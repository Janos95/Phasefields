//
// Created by janos on 16.12.19.
//

#pragma once

#include "IpTNLP.hpp"

#include <Eigen/Core>




template<class F, class... Derived>
void evaluateFunctional(const F& f, Ipopt::Number* x, const Ipopt::Index n, Eigen::PlainObjectBase<Derived>&... out)
{
    using VectorType = Eigen::Matrix<Ipopt::Number, Eigen::Dynamic, 1>;
    Eigen::Map<VectorType> mapped(x, n);
    f.evaluate(mapped, out...);
}

template<typename T, typename Enable = void>
class HessianBase
{
};


template<class F>
struct HessianBase<F, std::enable_if_t<provides_hessian_v<F>>>
{
    using HessianType = typename F::HessianType;
    HessianType h;
};



struct Solver
{
    struct Options{
        std::vector<std::function<void()>> callbacks;
    };

    struct Summary{

    };

    template<class P>
    Solver(const Options& options, P&& problem, Summary& summary);
};



template<class F, class G>
class IpoptProblem: public Ipopt::TNLP
{
    F& m_f;
    G& m_g;

    using DomainType = typename F::DomainType;
    using ResultType = typename F::ResultType;
    using GradientType = typename F::GradientType;

    ResultType m_fResult;
    GradientType m_fGradient;

    ResultType m_gResult;
    GradientType m_gGradient;

    HessianBase<F> m_fHessian;
    HessianBase<G> m_gHessian;

    DomainType m_x;

public:

    struct Options : Solver::Options{
        bool check_gradient = false;
    };

    IpoptProblem(F&& f, G&& g):
        m_f(std::forward<F>(f)),
        m_g(std::forward<G>(g))
    {
    }


private:


    /**@name Overloaded from TNLP */
    //@{
    /** Method to return some info about the NLP */
     bool get_nlp_info(
            Ipopt::Index&          n,
            Ipopt::Index&          m,
            Ipopt::Index&          nnz_jac_g,
            Ipopt::Index&          nnz_h_lag,
            IndexStyleEnum& index_style
    ) override
    {
        n = m_f.getDomainDimension();
        m = m_g.getCoDomainDimension();

        nnz_jac_g = m_g.getJacNNZ();
        nnz_jac_f = m_f.getJacNNZ();


        index_style = TNLP::C_STYLE;
    }

    /** Method to return the bounds for my problem */
     bool get_bounds_info(
            Ipopt::Index   n,
            Ipopt::Number* x_l,
            Ipopt::Number* x_u,
            Ipopt::Index   m,
            Ipopt::Number* g_l,
            Ipopt::Number* g_u
    ) override
    {

    }

    /** Method to return the starting point for the algorithm */
     bool get_starting_point(
            Ipopt::Index   n,
            bool    init_x,
            Ipopt::Number* x,
            bool    init_z,
            Ipopt::Number* z_L,
            Ipopt::Number* z_U,
            Ipopt::Index   m,
            bool    init_lambda,
            Ipopt::Number* lambda
    ) override;

    /** Method to return the objective value */
     bool eval_f(
            Ipopt::Index         n,
            const Ipopt::Number* x,
            bool          new_x,
            Ipopt::Number&       obj_value
    ) override
    {
         if(new_x){
             Eigen::Map<DomainType> mapped(x);
             if constexpr(provides_hessian_v<F>)
                 m_f.evaluate(mapped, m_fResult, m_fGradient, m_fHessian.h);
             else
                 m_f.evaluate(mapped, m_fResult, m_fGradient);
         }

         obj_value = m_fResult;
         return true;
    }

    /** Method to return the gradient of the objective */
     bool eval_grad_f(
            Ipopt::Index         n,
            const Ipopt::Number* x,
            bool          new_x,
            Ipopt::Number*       grad_f
    ) override
    {
        if(new_x){
            Eigen::Map<DomainType> mapped(x);
            if constexpr(provides_hessian_v<F>)
                m_f.evaluate(mapped, m_fResult, m_fGradient, m_fHessian.h);
            else
                m_f.evaluate(mapped, m_fResult, m_fGradient);
        }

        std::copy_n(m_fGradient.begin(), n, grad_f);
        return true;
    }


    /**
     * @brief no constraints for now so just return 1.
     */
     bool eval_g(
            Ipopt::Index         n,
            const Ipopt::Number* x,
            bool          new_x,
            Ipopt::Index         m,
            Ipopt::Number*       g
    ) override
    {
         return true;
    }

    /**
     * @brief no constraints for now so just return 1.
     */
    bool eval_jac_g(
            Ipopt::Index         n,
            const Ipopt::Number* x,
            bool          new_x,
            Ipopt::Index         m,
            Ipopt::Index         nele_jac,
            Ipopt::Index*        iRow,
            Ipopt::Index*        jCol,
            Ipopt::Number*       values
    ) override
    {
         return true;
    }

    /** Method to return:
     *   1) The structure of the hessian of the lagrangian (if "values" is NULL)
     *   2) The values of the hessian of the lagrangian (if "values" is not NULL)
     */
     bool eval_h(
            Ipopt::Index         n,
            const Ipopt::Number* x,
            bool          new_x,
            Ipopt::Number        obj_factor,
            Ipopt::Index         m,
            const Ipopt::Number* lambda,
            bool          new_lambda,
            Ipopt::Index         nele_hess,
            Ipopt::Index*        iRow,
            Ipopt::Index*        jCol,
            Ipopt::Number*       values
    ) override
    {
         return false;
    }

     void finalize_solution(
            Ipopt::SolverReturn               status,
            Ipopt::Index                      n,
            const Ipopt::Number*              x,
            const Ipopt::Number*              z_L,
            const Ipopt::Number*              z_U,
            Ipopt::Index                      m,
            const Ipopt::Number*              g,
            const Ipopt::Number*              lambda,
            Ipopt::Number                     obj_value,
            const Ipopt::IpoptData*           ip_data,
            Ipopt::IpoptCalculatedQuantities* ip_cq
    ) override
     {

     }





};

