//
// Created by janos on 20.04.20.
//

#include "Solver.h"
#include "RecursiveProblem.h"
#include "SparseMatrix.h"
#include "FunctionRef.h"

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
#include <Corrade/Containers/Pointer.h>

namespace Phasefield::Solver {

namespace {

struct IpoptWrapper : Ipopt::TNLP {
    const double Infinity = +2e19;
    const double NegativeInfinity = -2e19;

    explicit IpoptWrapper(
            Solver::RecursiveProblem const& pb,
            Solver::Options& opt,
            Containers::ArrayView<double> const& params) :
            problem(pb), options(opt), parameters(params),
            constraints(Containers::NoInit, problem.numConstraints()),
            gradient(Containers::NoInit, problem.numParameters()) {
    }

    Solver::RecursiveProblem const& problem;
    Solver::Options& options;

    Containers::ArrayView<double> parameters;
    double objective;
    Containers::Array<double> constraints;
    Containers::Array<double> gradient; /* gradient of the objective function */
    SparseMatrix jacobian;

    static_assert(std::is_same_v<Ipopt::Number, double>);

    void eval(double const* const x) {
        problem({x, problem.numParameters()}, objective, gradient, constraints, &jacobian);
    }

    bool get_nlp_info(
            Ipopt::Index& n,
            Ipopt::Index& m,
            Ipopt::Index& nnz_jac_g,
            Ipopt::Index& nnz_h_lag,
            IndexStyleEnum& index_style
    ) override {

        n = problem.numParameters();
        m = problem.numConstraints();
        problem.determineSparsityStructure(jacobian);

        nnz_jac_g = jacobian.nnz;

        /* 0 based indexing (i.e. how any sane human would index arrays */
        index_style = TNLP::C_STYLE;
        return true;
    }

    /** Method to return the bounds for my problem */
    bool get_bounds_info(
            Ipopt::Index n,
            Ipopt::Number* x_l,
            Ipopt::Number* x_u,
            Ipopt::Index m,
            Ipopt::Number* g_l,
            Ipopt::Number* g_u
    ) override {
        for(int i = 0; i < n; ++i){
            x_l[i] = NegativeInfinity;
            x_u[i] = Infinity;
        }
        for(int i = 0; i < m; ++i){
            g_l[i] = 0.;
            g_u[i] = 0.;
        }
        return true;
    }

    /** Method to return the starting point for the algorithm */
    bool get_starting_point(
            Ipopt::Index n,
            bool init_x,
            Ipopt::Number* x,
            bool init_z,
            Ipopt::Number* z_L,
            Ipopt::Number* z_U,
            Ipopt::Index m,
            bool init_lambda,
            Ipopt::Number* lambda
    ) override {
        // Here, we assume we only have starting values for x, if you code
        // your own NLP, you can provide starting values for the dual variables
        // if you wish
        CORRADE_INTERNAL_ASSERT(init_x == true);
        CORRADE_INTERNAL_ASSERT(init_z == false);
        CORRADE_INTERNAL_ASSERT(init_lambda == false);

        // initialize to the given starting point
        Cr::Utility::copy(parameters, {x, std::size_t(n)});
        return true;
    }

    /** Method to return the objective value */
    bool eval_f(
            Ipopt::Index n,
            const Ipopt::Number* x,
            bool new_x,
            Ipopt::Number& obj_value
    ) override {
        if(new_x) eval(x);
        obj_value = objective;
        return true;
    }

    /** Method to return the gradient of the objective */
    bool eval_grad_f(
            Ipopt::Index n,
            const Ipopt::Number* x,
            bool new_x,
            Ipopt::Number* grad_f
    ) override {
        if(new_x) eval(x);
        Cr::Utility::copy(gradient, {grad_f, (std::size_t) n});
        return true;
    }

    /**
     * @brief no constraints for now so just return 1.
     */
    bool eval_g(
            Ipopt::Index n,
            const Ipopt::Number* x,
            bool new_x,
            Ipopt::Index m,
            Ipopt::Number* g
    ) override {
        if(new_x) eval(x);
        Cr::Utility::copy(constraints, {g, (std::size_t) m});
        return true;
    }

    bool eval_jac_g(
            Ipopt::Index n,
            const Ipopt::Number* x,
            bool new_x,
            Ipopt::Index m,
            Ipopt::Index nele_jac,
            Ipopt::Index* iRow,
            Ipopt::Index* jCol,
            Ipopt::Number* values
    ) override {
        if(values){
            if(new_x) eval(x);
            Cr::Utility::copy(jacobian.values, {values, (std::size_t) nele_jac});
        } else{
            Cr::Utility::copy(jacobian.rows, {iRow, (std::size_t) nele_jac});
            Cr::Utility::copy(jacobian.cols, {jCol, (std::size_t) nele_jac});
        }
        return true;
    }

    /* only first order atm */
    bool eval_h(
            Ipopt::Index n,
            const Ipopt::Number* x,
            bool new_x,
            Ipopt::Number obj_factor,
            Ipopt::Index m,
            const Ipopt::Number* lambda,
            bool new_lambda,
            Ipopt::Index nele_hess,
            Ipopt::Index* iRow,
            Ipopt::Index* jCol,
            Ipopt::Number* values
    ) override {
        return false;
    }

    void finalize_solution(
            Ipopt::SolverReturn status,
            Ipopt::Index n,
            const Ipopt::Number* x,
            const Ipopt::Number* z_L,
            const Ipopt::Number* z_U,
            Ipopt::Index m,
            const Ipopt::Number* g,
            const Ipopt::Number* lambda,
            Ipopt::Number obj_value,
            const Ipopt::IpoptData* ip_data,
            Ipopt::IpoptCalculatedQuantities* ip_cq
    ) override {
        Cr::Utility::copy({x, (std::size_t) n}, parameters);
    }

    bool intermediate_callback(
            Ipopt::AlgorithmMode mode,
            Ipopt::Index iter,
            Ipopt::Number obj_value,
            Ipopt::Number inf_pr,
            Ipopt::Number inf_du,
            Ipopt::Number mu,
            Ipopt::Number d_norm,
            Ipopt::Number regularization_size,
            Ipopt::Number alpha_du,
            Ipopt::Number alpha_pr,
            Ipopt::Index ls_trials,
            const Ipopt::IpoptData* ip_data,
            Ipopt::IpoptCalculatedQuantities* ip_cq
    ) override {
        Solver::IterationSummary solverSummary; /* dummy variable */
        for(auto& cb : options.callbacks){
            switch(cb(solverSummary)) {
                case Solver::Status::CONTINUE :
                    return true;
                case Solver::Status::USER_ABORTED :
                case Solver::Status::FINISHED :
                    return false;
            }
        }
        return true;
    }
};

struct FirstOrderWrapper : ceres::FirstOrderFunction {

    Solver::RecursiveProblem& problem;
    mutable SparseMatrix jacobian;
    mutable Containers::Array<double> constraints;

    explicit FirstOrderWrapper(Solver::RecursiveProblem& pb) : problem(pb), constraints(Containers::NoInit, pb.numConstraints()) {
        CORRADE_ASSERT(pb.constraints.empty(), "Solver : ceres does not support constraints",);
        problem.determineSparsityStructure(jacobian);
    }


    bool Evaluate(double const* parameters, double* cost, double* g) const override {

        auto n = problem.numParameters();

        if(g)
            for(size_t i = 0; i < n; ++i) g[i] = 0.;
        *cost = 0;
        // @TODO where do we zero the data for sure?
        problem({parameters, n}, *cost, {g, g ? n : 0});

        if(g){
            for(Int i = 0; i < jacobian.nnz; ++i){
                g[jacobian.cols[i]] += jacobian.values[i];
            }
        }

        return true;
    }

    [[nodiscard]] int NumParameters() const override { return problem.numParameters(); }
};

struct CeresCallbackWrapper : ceres::IterationCallback {
    using callback_type = FunctionRef<Solver::Status::Value(Solver::IterationSummary const&)>;

    CeresCallbackWrapper(callback_type const& cb) : callback(cb) {}

    callback_type callback;

    ceres::CallbackReturnType operator()(const ceres::IterationSummary& summary) override {
        Solver::IterationSummary solverSummary;
        switch(callback(solverSummary)) {
            case Solver::Status::USER_ABORTED :
                return ceres::CallbackReturnType::SOLVER_ABORT;
            case Solver::Status::CONTINUE :
                return ceres::CallbackReturnType::SOLVER_CONTINUE;
            case Solver::Status::FINISHED :
                return ceres::CallbackReturnType::SOLVER_TERMINATE_SUCCESSFULLY;
        }
        CORRADE_ASSERT(false, "Optimization Callback : Unknown status", ceres::CallbackReturnType::SOLVER_ABORT);
    }
};

ceres::LineSearchDirectionType mapLineSearchType(Solver::LineSearchDirection::Value type) {
    switch(type) {
        case Solver::LineSearchDirection::NONLINEAR_CONJUGATE_GRADIENT :
            return ceres::NONLINEAR_CONJUGATE_GRADIENT;
        case Solver::LineSearchDirection::LBFGS :
            return ceres::LBFGS;
        case Solver::LineSearchDirection::BFGS :
            return ceres::BFGS;
        case Solver::LineSearchDirection::STEEPEST_DESCENT :
            return ceres::STEEPEST_DESCENT;
        default :
            CORRADE_ASSERT(false, "Line Search Type Not Supported", {});
    }
}

}

void solve(Solver::Options& options, Solver::RecursiveProblem& problem, Containers::ArrayView<double> params, Solver::Summary* summary) {
    if(options.solver == Solver::Backend::CERES) {
        ceres::GradientProblemSolver::Summary ceresSummary;
        ceres::GradientProblem ceresProblem(new FirstOrderWrapper(problem));
        ceres::GradientProblemSolver::Options ceresOptions{
                .line_search_direction_type = mapLineSearchType(options.line_search_direction),
                .max_num_iterations = int(options.max_num_iterations),
                .minimizer_progress_to_stdout = options.minimizer_progress_to_stdout,
                .update_state_every_iteration = options.update_state_every_iteration,
        };

        Array<Pointer<CeresCallbackWrapper>> cbs(options.callbacks.size());
        for(std::size_t i = 0; i < cbs.size(); ++i)
            cbs[i].reset(new CeresCallbackWrapper(options.callbacks[i]));

        for(auto& cb : cbs) ceresOptions.callbacks.push_back(cb.get());

        /* solve the problem and print the report */
        ceres::Solve(ceresOptions, ceresProblem, params.data(), &ceresSummary);
        Debug{} << ceresSummary.BriefReport().c_str();

    } else if(options.solver == Solver::Backend::IPOPT) {
        Ipopt::SmartPtr<Ipopt::TNLP> mynlp = new IpoptWrapper(problem, options, params);
        Ipopt::SmartPtr<Ipopt::IpoptApplication> app = IpoptApplicationFactory();
        app->Options()->SetNumericValue("tol", 1e-7);
        app->Options()->SetIntegerValue("print_level", 8);
        app->Options()->SetStringValue("hessian_approximation", "limited-memory");
        if(options.line_search_direction == Solver::LineSearchDirection::LBFGS)
            app->Options()->SetStringValue("limited_memory_update_type", "bfgs");
        else if(options.line_search_direction == Solver::LineSearchDirection::SR1)
            app->Options()->SetStringValue("limited_memory_update_type", "sr1");
        else
            CORRADE_ASSERT(false, "Line search not supported by ipopt",);
        Ipopt::ApplicationReturnStatus status = app->Initialize();
        if(status < 0) {
            Debug{} << "Ipopt : Error initializing solver";
            return;
        }

        status = app->OptimizeTNLP(mynlp);
        if(status < 0) {
            Debug{} << "Ipopt : Error solving problem";
        }

    } else
        CORRADE_ASSERT(false, "Unkown solver type",);
}

}


