//
// Created by janos on 19.11.19.
//

#pragma once

#include <IpTNLP.hpp>

#include <Eigen/Core>

class Optimization : public Ipopt::TNLP
{

private:

    Eigen::VectorXd phaseField;

    bool get_nlp_info(Ipopt::Index& n, Ipopt::Index& m, Ipopt::Index& nnz_jac_g,
                              Ipopt::Index& nnz_h_lag, IndexStyleEnum& index_style) override
    {
        n = phaseField.cols();
        m = 1;

        nnz_jac_g = -1; //TODO
        nnz_h_lag = -1; //TODO


        index_style = TNLP::C_STYLE;

        return true;
    }

    bool get_bounds_info(
            Ipopt::Index n, Ipopt::Number* x_l, Ipopt::Number* x_u,
            Ipopt::Index m, Ipopt::Number* g_l, Ipopt::Number* g_u) override
    {

        // the variables have lower bounds of 1
        std::fill_n(x_l, x_l + n, Ipopt::nlp_lower_bound_inf);
        std::fill_n(x_u, x_u + n, Ipopt::nlp_upper_bound_inf);


        return true;
    }


};
