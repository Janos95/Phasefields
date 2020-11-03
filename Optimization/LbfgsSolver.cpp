//
// Created by janos on 11/1/20.
//

#include "LbfgsSolver.h"
#include "Solver.h"

#include <limits>

namespace Phasefield {

/* yeah this is a bit hacky .. */
namespace {
#include "../contrib/liblbfgs/lib/lbfgs.c"


struct lbfgs_data {
    Solver::RecursiveProblem& problem;
    Solver::Options& options;
};


lbfgsfloatval_t evaluate(
        void *instance,
        const lbfgsfloatval_t *x,
        lbfgsfloatval_t *g,
        const int n,
        const lbfgsfloatval_t step
)
{
    size_t N = n;
    auto& problem = reinterpret_cast<lbfgs_data*>(instance)->problem;
    double fx = 0;
    problem({x, N}, fx, {g, N});
    return fx;
}


int checkInput(lbfgs_parameter_t& param, line_search_proc& linesearch, size_t n) {
    /* Check the input parameters for errors. */
    if (n <= 0) {
        return LBFGSERR_INVALID_N;
    }

    if (param.epsilon < 0.) {
        return LBFGSERR_INVALID_EPSILON;
    }
    if (param.past < 0) {
        return LBFGSERR_INVALID_TESTPERIOD;
    }
    if (param.delta < 0.) {
        return LBFGSERR_INVALID_DELTA;
    }
    if (param.min_step < 0.) {
        return LBFGSERR_INVALID_MINSTEP;
    }
    if (param.max_step < param.min_step) {
        return LBFGSERR_INVALID_MAXSTEP;
    }
    if (param.ftol < 0.) {
        return LBFGSERR_INVALID_FTOL;
    }
    if (param.linesearch == LBFGS_LINESEARCH_BACKTRACKING_WOLFE ||
        param.linesearch == LBFGS_LINESEARCH_BACKTRACKING_STRONG_WOLFE) {
        if (param.wolfe <= param.ftol || 1. <= param.wolfe) {
            return LBFGSERR_INVALID_WOLFE;
        }
    }
    if (param.gtol < 0.) {
        return LBFGSERR_INVALID_GTOL;
    }
    if (param.xtol < 0.) {
        return LBFGSERR_INVALID_XTOL;
    }
    if (param.max_linesearch <= 0) {
        return LBFGSERR_INVALID_MAXLINESEARCH;
    }
    if (param.orthantwise_c < 0.) {
        return LBFGSERR_INVALID_ORTHANTWISE;
    }
    if (param.orthantwise_start < 0 || n < param.orthantwise_start) {
        return LBFGSERR_INVALID_ORTHANTWISE_START;
    }
    if (param.orthantwise_end < 0) {
        param.orthantwise_end = n;
    }
    if (n < param.orthantwise_end) {
        return LBFGSERR_INVALID_ORTHANTWISE_END;
    }
    if (param.orthantwise_c != 0.) {
        switch (param.linesearch) {
            case LBFGS_LINESEARCH_BACKTRACKING:
                linesearch = line_search_backtracking_owlqn;
                break;
            default:
                /* Only the backtracking method is available. */
                return LBFGSERR_INVALID_LINESEARCH;
        }
    } else {
        switch (param.linesearch) {
            case LBFGS_LINESEARCH_MORETHUENTE:
                linesearch = line_search_morethuente;
                break;
            case LBFGS_LINESEARCH_BACKTRACKING_ARMIJO:
            case LBFGS_LINESEARCH_BACKTRACKING_WOLFE:
            case LBFGS_LINESEARCH_BACKTRACKING_STRONG_WOLFE:
                linesearch = line_search_backtracking;
                break;
            default:
                return LBFGSERR_INVALID_LINESEARCH;
        }
    }

    return -1;
}

static int progress(
        void *instance,
        const lbfgsfloatval_t *x,
        const lbfgsfloatval_t *g,
        const lbfgsfloatval_t fx,
        const lbfgsfloatval_t xnorm,
        const lbfgsfloatval_t gnorm,
        const lbfgsfloatval_t step,
        int n,
        int k,
        int ls
)
{
    auto& options = reinterpret_cast<lbfgs_data*>(instance)->options;
    Solver::IterationSummary solverSummary; /* dummy variable */
    bool goOn = true;
    for(auto& cb : options.callbacks){
        switch(cb(solverSummary)) {
            case Solver::Status::CONTINUE :
                break;
            case Solver::Status::USER_ABORTED :
            case Solver::Status::FINISHED :
                goOn = false;
        }
    }
    return goOn ? 0 : 1;
}

}


struct LbfgsSolver::Data {
    int ret;
    int i, j, k, ls, end, bound;
    lbfgsfloatval_t step;

    /* Constant parameters and their default values. */
    lbfgs_parameter_t param;
    int m;

    lbfgsfloatval_t *xp = NULL;
    lbfgsfloatval_t *g = NULL, *gp = NULL, *pg = NULL;
    lbfgsfloatval_t *d = NULL, *w = NULL, *pf = NULL;
    iteration_data_t *lm = NULL, *it = NULL;
    lbfgsfloatval_t ys, yy;
    lbfgsfloatval_t xnorm, gnorm, beta;
    lbfgsfloatval_t fx = 0.;
    lbfgsfloatval_t rate = 0.;
    line_search_proc linesearch = line_search_morethuente;
};

LbfgsSolver::LbfgsSolver(Solver::Options& options, Solver::RecursiveProblem& problem, ArrayView<double> data) :
    m_data(new LbfgsSolver::Data{}),
    m_cd(new callback_data_t{}),
    m_problem(problem),
    m_options(options),
    m_x(data)
{
    lbfgs_parameter_t param;

    /* Initialize the parameters for the L-BFGS optimization. */
    lbfgs_parameter_init(&param);
    param.epsilon = 1e-10;
    param.ftol = 1e-10;
    param.max_iterations = options.max_num_iterations;
    param.xtol = std::numeric_limits<double>::epsilon();
    param.delta = 1e-10;

    size_t m = m_data->m;
    size_t n = data.size();

    /* Allocate working space. */
    m_data->xp = (lbfgsfloatval_t*)vecalloc(n * sizeof(lbfgsfloatval_t));
    m_data->g = (lbfgsfloatval_t*)vecalloc(n * sizeof(lbfgsfloatval_t));
    m_data->gp = (lbfgsfloatval_t*)vecalloc(n * sizeof(lbfgsfloatval_t));
    m_data->d = (lbfgsfloatval_t*)vecalloc(n * sizeof(lbfgsfloatval_t));
    m_data->w = (lbfgsfloatval_t*)vecalloc(n * sizeof(lbfgsfloatval_t));

    if (param.orthantwise_c != 0.) {
        /* Allocate working space for OW-LQN. */
        m_data->pg = (lbfgsfloatval_t*)vecalloc(n * sizeof(lbfgsfloatval_t));
    }

    /* Allocate limited memory storage. */
    m_data->lm = (iteration_data_t*)vecalloc(m * sizeof(iteration_data_t));

    /* Initialize the limited memory. */
    for (size_t i = 0; i < m;++i) {
        auto it = &(m_data->lm[i]);
        it->alpha = 0;
        it->ys = 0;
        it->s = (lbfgsfloatval_t*)vecalloc(n * sizeof(lbfgsfloatval_t));
        it->y = (lbfgsfloatval_t*)vecalloc(n * sizeof(lbfgsfloatval_t));
    }

    /* Allocate an array for storing previous values of the objective function. */
    if (0 < param.past) {
        m_data->pf = (lbfgsfloatval_t*)vecalloc(param.past * sizeof(lbfgsfloatval_t));
    }

    checkInput(m_data->param, m_data->linesearch, n);

    auto* x = m_x.data();
    auto g = m_data->g;
    auto pg = m_data->pg;
    auto pf = m_data->pf;
    auto d = m_data->d;

    /* Construct a callback data. */
    auto& cd = *(callback_data_t*)m_cd;
    lbfgs_data* instance = new lbfgs_data{m_problem, m_options};
    cd.n = n;
    cd.instance = instance;
    cd.proc_evaluate = evaluate;
    cd.proc_progress = progress;

    /* Evaluate the function value and its gradient. */
    fx = cd.proc_evaluate(cd.instance, x, g, cd.n, 0);
    if (0. != param.orthantwise_c) {
        /* Compute the L1 norm of the variable and add it to the object value. */
        xnorm = owlqn_x1norm(x, param.orthantwise_start, param.orthantwise_end);
        fx += xnorm * param.orthantwise_c;
        owlqn_pseudo_gradient(
                pg, x, g, n,
                param.orthantwise_c, param.orthantwise_start, param.orthantwise_end
        );
    }

    /* Store the initial value of the objective function. */
    if (pf != NULL) {
        pf[0] = fx;
    }

    /*
        Compute the direction;
        we assume the initial hessian matrix H_0 as the identity matrix.
     */
    if (param.orthantwise_c == 0.) {
        vecncpy(d, g, n);
    } else {
        vecncpy(d, pg, n);
    }

    /*
       Make sure that the initial variables are not a minimizer.
     */
    vec2norm(&xnorm, x, n);
    if (param.orthantwise_c == 0.) {
        vec2norm(&gnorm, g, n);
    } else {
        vec2norm(&gnorm, pg, n);
    }
    if (xnorm < 1.0) xnorm = 1.0;
    if (gnorm / xnorm <= param.epsilon) {
        ret = LBFGS_ALREADY_MINIMIZED;
    }

    /* Compute the initial step:
        step = 1.0 / sqrt(vecdot(d, d, n))
     */
    vec2norminv(&step, d, n);

    k = 1;
    end = 0;



}

LbfgsSolver::~LbfgsSolver() {
    vecfree(m_data->pf);

    /* Free memory blocks used by this function. */
    if (m_data->lm != NULL) {
        auto lm = m_data->lm;
        for (int i = 0;i < m_data->m;++i) {
            vecfree(lm[i].s);
            vecfree(lm[i].y);
        }
        vecfree(lm);
    }
    vecfree(m_data->pg);
    vecfree(m_data->w);
    vecfree(m_data->d);
    vecfree(m_data->gp);
    vecfree(m_data->g);
    vecfree(m_data->xp);

    delete m_data;
    auto* cd = (callback_data_t*)m_cd;
    delete cd->instance;
    delete cd;
}


int LbfgsSolver::runOneIteration() {
    auto& cd = *(callback_data_t*)m_cd;
    lbfgs_parameter_t param = m_data->param;

    line_search_proc linesearch = m_data->linesearch;

    lbfgsfloatval_t* x = m_x.data();

    const size_t m = param.m;
    const size_t n = m_x.size();

    lbfgsfloatval_t *xp = m_data->xp;
    lbfgsfloatval_t *g = m_data->g, *gp = m_data->gp, *pg = m_data->pg;
    lbfgsfloatval_t *d = m_data->d, *w = m_data->w, *pf = m_data->pf;
    iteration_data_t *lm = m_data->lm, *it = m_data->it;

    /* ------------------------------------------------------- */


    //for (;;) {
        /* Store the current position and gradient vectors. */
        veccpy(xp, x, n);
        veccpy(gp, g, n);

        /* Search for an optimal step. */
        if (param.orthantwise_c == 0.) {
            ls = linesearch(n, x, &fx, g, d, &step, xp, gp, w, &cd, &param);
        } else {
            ls = linesearch(n, x, &fx, g, d, &step, xp, pg, w, &cd, &param);
            owlqn_pseudo_gradient(
                    pg, x, g, n,
                    param.orthantwise_c, param.orthantwise_start, param.orthantwise_end
            );
        }
        if (ls < 0) {
            /* Revert to the previous point. */
            veccpy(x, xp, n);
            veccpy(g, gp, n);
            ret = ls;
            return ret;
        }

        /* Compute x and g norms. */
        vec2norm(&xnorm, x, n);
        if (param.orthantwise_c == 0.) {
            vec2norm(&gnorm, g, n);
        } else {
            vec2norm(&gnorm, pg, n);
        }

        /* Report the progress. */
        if (cd.proc_progress) {
            if ((ret = cd.proc_progress(cd.instance, x, g, fx, xnorm, gnorm, step, cd.n, k, ls))) {
                return ret;
            }
        }

        /*
            Convergence test.
            The criterion is given by the following formula:
                |g(x)| / \max(1, |x|) < \epsilon
         */
        if (xnorm < 1.0) xnorm = 1.0;
        if (gnorm / xnorm <= param.epsilon) {
            /* Convergence. */
            return LBFGS_SUCCESS;
        }

        /*
            Test for stopping criterion.
            The criterion is given by the following formula:
                |(f(past_x) - f(x))| / f(x) < \delta
         */
        if (pf != NULL) {
            /* We don't test the stopping criterion while k < past. */
            if (param.past <= k) {
                /* Compute the relative improvement from the past. */
                rate = (pf[k % param.past] - fx) / fx;

                /* The stopping criterion. */
                if (fabs(rate) < param.delta) {
                    return LBFGS_STOP;
                }
            }

            /* Store the current value of the objective function. */
            pf[k % param.past] = fx;
        }

        if (param.max_iterations != 0 && param.max_iterations < k+1) {
            /* Maximum number of iterations. */
            return LBFGSERR_MAXIMUMITERATION;
        }

        /*
            Update vectors s and y:
                s_{k+1} = x_{k+1} - x_{k} = \step * d_{k}.
                y_{k+1} = g_{k+1} - g_{k}.
         */
        it = &lm[end];
        vecdiff(it->s, x, xp, n);
        vecdiff(it->y, g, gp, n);

        /*
            Compute scalars ys and yy:
                ys = y^t \cdot s = 1 / \rho.
                yy = y^t \cdot y.
            Notice that yy is used for scaling the hessian matrix H_0 (Cholesky factor).
         */
        vecdot(&ys, it->y, it->s, n);
        vecdot(&yy, it->y, it->y, n);
        it->ys = ys;

        /*
            Recursive formula to compute dir = -(H \cdot g).
                This is described in page 779 of:
                Jorge Nocedal.
                Updating Quasi-Newton Matrices with Limited Storage.
                Mathematics of Computation, Vol. 35, No. 151,
                pp. 773--782, 1980.
         */
        bound = (m <= k) ? m : k;
        ++k;
        end = (end + 1) % m;

        /* Compute the steepest direction. */
        if (param.orthantwise_c == 0.) {
            /* Compute the negative of gradients. */
            vecncpy(d, g, n);
        } else {
            vecncpy(d, pg, n);
        }

        j = end;
        for (i = 0;i < bound;++i) {
            j = (j + m - 1) % m;    /* if (--j == -1) j = m-1; */
            it = &lm[j];
            /* \alpha_{j} = \rho_{j} s^{t}_{j} \cdot q_{k+1}. */
            vecdot(&it->alpha, it->s, d, n);
            it->alpha /= it->ys;
            /* q_{i} = q_{i+1} - \alpha_{i} y_{i}. */
            vecadd(d, it->y, -it->alpha, n);
        }

        vecscale(d, ys / yy, n);

        for (i = 0;i < bound;++i) {
            it = &lm[j];
            /* \beta_{j} = \rho_{j} y^t_{j} \cdot \gamma_{i}. */
            vecdot(&beta, it->y, d, n);
            beta /= it->ys;
            /* \gamma_{i+1} = \gamma_{i} + (\alpha_{j} - \beta_{j}) s_{j}. */
            vecadd(d, it->s, it->alpha - beta, n);
            j = (j + 1) % m;        /* if (++j == m) j = 0; */
        }

        /*
            Constrain the search direction for orthant-wise updates.
         */
        if (param.orthantwise_c != 0.) {
            for (i = param.orthantwise_start;i < param.orthantwise_end;++i) {
                if (d[i] * pg[i] >= 0) {
                    d[i] = 0;
                }
            }
        }

        /*
            Now the search direction d is ready. We try step = 1 first.
         */
        step = 1.0;

        return ret;
    //}
}



}