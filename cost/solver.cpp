//
// Created by janos on 20.04.20.
//

#include "solver.hpp"

#include <ceres/iteration_callback.h>
#include <ceres/gradient_problem_solver.h>
#include <ceres/gradient_problem.h>
#include <ceres/first_order_function.h>

#include "IpTNLP.hpp"

using namespace Corrade;

struct IpoptWrapper;

struct FirstOrderWrapper : ceres::FirstOrderFunction {
    explicit FirstOrderWrapper(solver::Problem const& pb) : problem(pb) {
        CORRADE_ASSERT(pb.constraints.empty(), "Solver : ceres does not support constraints",);
    }

    solver::Problem const& problem;

    bool Evaluate(double const* parameters, double* cost, double* jacobian) const override {
        return problem.evaluate(parameters, cost, jacobian);
    }

    [[nodiscard]] int NumParameters() const override { return problem.numParameters(); }
};

struct CeresCallbackWrapper : ceres::IterationCallback {
    using callback_type = unique_function<solver::Status(solver::IterationSummary const&)>;

    CeresCallbackWrapper(callback_type& cb) : callback(cb) {}

    callback_type& callback;
    ceres::CallbackReturnType operator()(const ceres::IterationSummary& summary) override {
        solver::IterationSummary solverSummary;
        switch (callback(solverSummary)) {
            case solver::Status::ABORT : return ceres::CallbackReturnType::SOLVER_ABORT;
            case solver::Status::CONTINUE : return ceres::CallbackReturnType::SOLVER_CONTINUE;
            case solver::Status::FINISHED : return ceres::CallbackReturnType::SOLVER_TERMINATE_SUCCESSFULLY;
        }
    }
};

void solve(solver::Options& options, solver::Problem& problem, double* params, solver::Summary* summary){
    if(options.solver == solver::Solver::CERES){
        ceres::GradientProblemSolver::Summary ceresSummary;
        ceres::GradientProblem ceresProblem(new FirstOrderWrapper(problem));
        ceres::GradientProblemSolver::Options ceresOptions{
            .max_num_iterations = options.max_num_iterations,
            .minimizer_progress_to_stdout = options.minimizer_progress_to_stdout
        };
        Containers::Array<Containers::Pointer<CeresCallbackWrapper>> cbs(options.callbacks.size());
        for (int i = 0; i < cbs.size(); ++i)
            cbs[i].reset(new CeresCallbackWrapper(options.callbacks[i]));
        for(auto& cb : cbs) ceresOptions.callbacks.push_back(cb.get());
        ceres::Solve(ceresOptions, ceresProblem, params, &ceresSummary);
    } else if(options.solver == solver::Solver::IPOPT){
        IpoptWrapper i
    } else CORRADE_ASSERT(false, "Unkown solver type", );
}



struct IpoptWrapper : Ipopt::TNLP
{

    explicit IpoptWrapper(solver::Problem const& pb) : problem(pb) {}

    solver::Problem const& problem;

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
    ) override {

    }

    /** Method to return the objective value */
    bool eval_f(
            Ipopt::Index         n,
            const Ipopt::Number* x,
            bool          new_x,
            Ipopt::Number&       obj_value
    ) override
    {

    }

    /** Method to return the gradient of the objective */
    bool eval_grad_f(
            Ipopt::Index         n,
            const Ipopt::Number* x,
            bool          new_x,
            Ipopt::Number*       grad_f
    ) override
    {

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
