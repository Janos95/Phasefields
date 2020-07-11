//
// Created by janos on 20.04.20.
//

#include "solver.hpp"
#include "problem.hpp"

#include <ceres/iteration_callback.h>
#include <ceres/gradient_problem_solver.h>
#include <ceres/gradient_problem.h>
#include <ceres/first_order_function.h>


#include <IpTNLP.hpp>
#include <IpIpoptCalculatedQuantities.hpp>
#include <IpIpoptData.hpp>
#include <IpTNLPAdapter.hpp>
#include <IpOrigIpoptNLP.hpp>
#include <IpIpoptApplication.hpp>

#include <Corrade/Utility/Debug.h>
#include <Corrade/Utility/Algorithms.h>
#include <Corrade/Containers/Array.h>


using namespace Corrade;

namespace {

    struct IpoptWrapper : Ipopt::TNLP
    {
        const double Infinity = +2e19;
        const double NegativeInfinity = -2e19;

        explicit IpoptWrapper(
                solver::Problem const& pb,
                solver::Options const& opt,
                Containers::ArrayView<double> const& params) :
                    problem(pb), options(opt), parameters(params),
                    constraints(Containers::NoInit, problem.numConstraints()),
                    gradientConstraints(Containers::NoInit, problem.numConstraints() * problem.numParameters()),
                    gradientObjective(Containers::NoInit, problem.numParameters())
        {
        }

        solver::Problem const& problem;
        solver::Options const& options;

        Containers::ArrayView<double> parameters;

        Containers::Array<double> constraints;
        Containers::Array<double> gradientConstraints;

        double objective;
        Containers::Array<double> gradientObjective;

        static_assert(std::is_same_v<Ipopt::Number, double>);

        void eval(double const* const x){
            problem.evaluate(x, &objective, gradientObjective.data(), constraints.data(), gradientConstraints.data());
        }

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

            n = problem.numParameters();
            m = problem.constraints.size();

            /* dense jacobian for the beginning */
            nnz_jac_g = problem.numParameters() * problem.numConstraints();

            nnz_h_lag = 0; //not using hessian atm

            // use the C style indexing (0-based)
            index_style = TNLP::C_STYLE;

            return true;
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
            for (int i = 0; i < n; ++i) {
                x_l[i] = NegativeInfinity;
                x_u[i] = Infinity;
            }
            for (int i = 0; i < m; ++i) {
                g_l[i] = 0.;
                g_u[i] = 0.;
            }
            return true;
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
            // Here, we assume we only have starting values for x, if you code
            // your own NLP, you can provide starting values for the dual variables
            // if you wish
            CORRADE_INTERNAL_ASSERT(init_x == true);
            CORRADE_INTERNAL_ASSERT(init_z == false);
            CORRADE_INTERNAL_ASSERT(init_lambda == false);

            // initialize to the given starting point
            Utility::copy(parameters, {x,(std::size_t)n});
            return true;
        }

        /** Method to return the objective value */
        bool eval_f(
                Ipopt::Index         n,
                const Ipopt::Number* x,
                bool          new_x,
                Ipopt::Number&       obj_value
        ) override
        {
            if(new_x) eval(x);
            obj_value = objective;
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
            if(new_x) eval(x);
            Utility::copy(gradientObjective, {grad_f, (std::size_t)n});
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
            if(new_x) eval(x);
            Utility::copy(constraints, {g, (std::size_t)m});
            return true;
        }

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
            CORRADE_INTERNAL_ASSERT(n*m == nele_jac);
            if(values){
                if(new_x) eval(x);
                Utility::copy(gradientConstraints, {values, (std::size_t)nele_jac});
            } else {
                /* our thing is dense */
                for (int i = 0; i < n; ++i) {
                    for (int j = 0; j < m; ++j) {
                        iRow[j * n + i] = i;
                        jCol[j * n + i] = j;
                    }
                }
            }
            return true;
        }

        /* only first order atm */
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
            Utility::copy({x,(std::size_t)n}, parameters);
        }

        bool intermediate_callback(
                Ipopt::AlgorithmMode              mode,
                Ipopt::Index                      iter,
                Ipopt::Number                     obj_value,
                Ipopt::Number                     inf_pr,
                Ipopt::Number                     inf_du,
                Ipopt::Number                     mu,
                Ipopt::Number                     d_norm,
                Ipopt::Number                     regularization_size,
                Ipopt::Number                     alpha_du,
                Ipopt::Number                     alpha_pr,
                Ipopt::Index                      ls_trials,
                const Ipopt::IpoptData*           ip_data,
                Ipopt::IpoptCalculatedQuantities* ip_cq
        ) override
        {
            solver::IterationSummary solverSummary; /* dummy variable */
            for (auto const& cb : options.callbacks) {
                switch (cb(solverSummary)) {
                    case solver::Status::CONTINUE : return true;
                    case solver::Status::ABORT :
                    case solver::Status::FINISHED : return false;
                }
            }
            return true;
        }
    };

    struct FirstOrderWrapper : ceres::FirstOrderFunction {
        explicit FirstOrderWrapper(solver::Problem const& pb) : problem(pb) {
            CORRADE_ASSERT(pb.constraints.empty(), "Solver : ceres does not support constraints",);
        }

        solver::Problem const& problem;

        bool Evaluate(double const* parameters, double* cost, double* jacobian) const override {
            return problem.evaluate(parameters, cost, jacobian, nullptr, nullptr);
        }

        [[nodiscard]] int NumParameters() const override { return problem.numParameters(); }
    };

    struct CeresCallbackWrapper : ceres::IterationCallback {
        using callback_type = function_ref<solver::Status::Value(solver::IterationSummary const&)>;

        CeresCallbackWrapper(callback_type& cb) : callback(cb) {}

        callback_type callback;
        ceres::CallbackReturnType operator()(const ceres::IterationSummary& summary) override {
            solver::IterationSummary solverSummary;
            switch (callback(solverSummary)) {
                case solver::Status::ABORT : return ceres::CallbackReturnType::SOLVER_ABORT;
                case solver::Status::CONTINUE : return ceres::CallbackReturnType::SOLVER_CONTINUE;
                case solver::Status::FINISHED : return ceres::CallbackReturnType::SOLVER_TERMINATE_SUCCESSFULLY;
            }
            CORRADE_ASSERT(false, "Optimization Callback : Unknown status", ceres::CallbackReturnType::SOLVER_ABORT);
        }
    };

    ceres::LineSearchDirectionType mapLineSearchType(solver::LineSearchDirection::Value type){
        switch (type) {
            case solver::LineSearchDirection::NONLINEAR_CONJUGATE_GRADIENT : return ceres::NONLINEAR_CONJUGATE_GRADIENT;
            case solver::LineSearchDirection::LBFGS : return ceres::LBFGS;
            case solver::LineSearchDirection::BFGS : return ceres::BFGS;
            case solver::LineSearchDirection::STEEPEST_DESCENT : return ceres::STEEPEST_DESCENT;
            default : CORRADE_ASSERT(false, "Line Search Type Not Supported", {});
        }
    }
}


void solve(solver::Options& options, solver::Problem& problem, double* params, solver::Summary* summary){
    if(options.solver == solver::Solver::CERES){
        ceres::GradientProblemSolver::Summary ceresSummary;
        ceres::GradientProblem ceresProblem(new FirstOrderWrapper(problem));
        ceres::GradientProblemSolver::Options ceresOptions{
            .line_search_direction_type = mapLineSearchType(options.line_search_direction),
            .max_num_iterations = options.max_num_iterations,
            .minimizer_progress_to_stdout = options.minimizer_progress_to_stdout,
            .update_state_every_iteration = options.update_state_every_iteration,
        };
        Containers::Array<Containers::Pointer<CeresCallbackWrapper>> cbs(options.callbacks.size());
        for (int i = 0; i < cbs.size(); ++i)
            cbs[i].reset(new CeresCallbackWrapper(options.callbacks[i]));
        for(auto& cb : cbs) ceresOptions.callbacks.push_back(cb.get());
        ceres::Solve(ceresOptions, ceresProblem, params, &ceresSummary);
        Mg::Debug{} << ceresSummary.BriefReport().c_str();
    } else if(options.solver == solver::Solver::IPOPT){
        Ipopt::SmartPtr<Ipopt::TNLP> mynlp = new IpoptWrapper(problem, options, {params, (std::size_t)problem.numParameters()});
        Ipopt::SmartPtr<Ipopt::IpoptApplication> app = IpoptApplicationFactory();
        app->Options()->SetNumericValue("tol", 1e-7);
        app->Options()->SetIntegerValue("print_level", 8);
        app->Options()->SetStringValue("hessian_approximation", "limited-memory");
        if(options.line_search_direction == solver::LineSearchDirection::LBFGS)
            app->Options()->SetStringValue("limited_memory_update_type", "bfgs");
        else if(options.line_search_direction == solver::LineSearchDirection::SR1)
            app->Options()->SetStringValue("limited_memory_update_type", "sr1");
        else CORRADE_ASSERT(false, "Line search not supported by ipopt",);
        auto status = app->Initialize();
        if( status != Ipopt::Solve_Succeeded )
            return;

        status = app->OptimizeTNLP(mynlp);

    } else CORRADE_ASSERT(false, "Unkown solver type",);
}




