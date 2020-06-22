//
// Created by janos on 04.06.20.
//

#pragma once

//quadrature rules taken from https://people.sc.fsu.edu/~jburkardt/datasets/quadrature_rules_tri/quadrature_rules_tri.html
//presumably they took them from Gilbert Strang, George Fix,
//An Analysis of the Finite Element Method

namespace {

    struct Point {
        double x, y;
    };

    template<int Order>
    constexpr Point Abscissas[Order] = {};

    template<int Order>
    constexpr double Weights[Order] = {};

    template<>
    constexpr Point Abscissas<3>[3] = {
            Point{0.50000000000000000000, 0.00000000000000000000},
            Point{0.50000000000000000000, 0.50000000000000000000},
            Point{0.00000000000000000000, 0.50000000000000000000}
    };

    template<>
    constexpr double Weights<3>[3] = {
            0.33333333333333333333,
            0.33333333333333333333,
            0.33333333333333333333
    };

    template<>
    constexpr Point Abscissas<4>[4] = {
            Point{0.33333333333333333333, 0.33333333333333333333},
            Point{0.60000000000000000000, 0.20000000000000000000},
            Point{0.20000000000000000000, 0.60000000000000000000},
            Point{0.20000000000000000000, 0.20000000000000000000},
    };

    template<>
    constexpr double Weights<4>[4] = {
            -0.56250000000000000000,
            0.52083333333333333333,
            0.52083333333333333333,
            0.52083333333333333333
    };

    template<>
    constexpr Point Abscissas<6>[6] = {
            Point{0.816847572980459, 0.091576213509771},
            Point{0.091576213509771, 0.816847572980459},
            Point{0.091576213509771, 0.091576213509771},
            Point{0.108103018168070, 0.445948490915965},
            Point{0.445948490915965, 0.108103018168070},
            Point{0.445948490915965, 0.445948490915965}
    };

    template<>
    constexpr double Weights<6>[6] = {
            0.109951743655322,
            0.109951743655322,
            0.109951743655322,
            0.223381589678011,
            0.223381589678011,
            0.223381589678011
    };

    template<>
    constexpr Point Abscissas<7>[7] = {
            Point{0.33333333333333333, 0.33333333333333333},
            Point{0.79742698535308720, 0.10128650732345633},
            Point{0.10128650732345633, 0.79742698535308720},
            Point{0.10128650732345633, 0.10128650732345633},
            Point{0.05971587178976981, 0.47014206410511505},
            Point{0.47014206410511505, 0.05971587178976981},
            Point{0.47014206410511505, 0.47014206410511505}
    };

    template<>
    constexpr double Weights<7>[7] = {
            0.22500000000000000,
            0.12593918054482717,
            0.12593918054482717,
            0.12593918054482717,
            0.13239415278850616,
            0.13239415278850616,
            0.13239415278850616
    };

    template<>
    constexpr Point Abscissas<9>[9] = {
            Point{0.124949503233232, 0.437525248383384},
            Point{0.437525248383384, 0.124949503233232},
            Point{0.437525248383384, 0.437525248383384},
            Point{0.797112651860071, 0.165409927389841},
            Point{0.797112651860071, 0.037477420750088},
            Point{0.165409927389841, 0.797112651860071},
            Point{0.165409927389841, 0.037477420750088},
            Point{0.037477420750088, 0.797112651860071},
            Point{0.037477420750088, 0.165409927389841}
    };

    template<>
    constexpr double Weights<9>[9] = {
            0.205950504760887,
            0.205950504760887,
            0.205950504760887,
            0.063691414286223,
            0.063691414286223,
            0.063691414286223,
            0.063691414286223,
            0.063691414286223,
            0.063691414286223
    };

    template<>
    constexpr Point Abscissas<12>[12] = {
            Point{0.873821971016996, 0.063089014491502},
            Point{0.063089014491502, 0.873821971016996},
            Point{0.063089014491502, 0.063089014491502},
            Point{0.501426509658179, 0.249286745170910},
            Point{0.249286745170910, 0.501426509658179},
            Point{0.249286745170910, 0.249286745170910},
            Point{0.636502499121399, 0.310352451033785},
            Point{0.636502499121399, 0.053145049844816},
            Point{0.310352451033785, 0.636502499121399},
            Point{0.310352451033785, 0.053145049844816},
            Point{0.053145049844816, 0.636502499121399},
            Point{0.053145049844816, 0.310352451033785}
    };

    template<>
    constexpr double Weights<12>[12] = {
            0.050844906370207,
            0.050844906370207,
            0.050844906370207,
            0.116786275726379,
            0.116786275726379,
            0.116786275726379,
            0.082851075618374,
            0.082851075618374,
            0.082851075618374,
            0.082851075618374,
            0.082851075618374,
            0.082851075618374
    };

    template<>
    constexpr Point Abscissas<13>[13] = {
            Point{0.333333333333333, 0.333333333333333},
            Point{0.479308067841923, 0.260345966079038},
            Point{0.260345966079038, 0.479308067841923},
            Point{0.260345966079038, 0.260345966079038},
            Point{0.869739794195568, 0.065130102902216},
            Point{0.065130102902216, 0.869739794195568},
            Point{0.065130102902216, 0.065130102902216},
            Point{0.638444188569809, 0.312865496004875},
            Point{0.638444188569809, 0.048690315425316},
            Point{0.312865496004875, 0.638444188569809},
            Point{0.312865496004875, 0.048690315425316},
            Point{0.048690315425316, 0.638444188569809},
            Point{0.048690315425316, 0.312865496004875}
    };

    template<>
    constexpr double Weights<13>[13] = {
            -0.149570044467670,
            0.175615257433204,
            0.175615257433204,
            0.175615257433204,
            0.053347235608839,
            0.053347235608839,
            0.053347235608839,
            0.077113760890257,
            0.077113760890257,
            0.077113760890257,
            0.077113760890257,
            0.077113760890257,
            0.077113760890257
    };
}
/**
 * integrates a function of type F over the
 * reference triangle {(0,0), (1,0), (0,1)}.
 * The function f is expected to take barycentric coordinates,
 * thus the integral of f over an arbitrary triangle T is 2 * |T| * integrate(f).
 * The factor of 2 stems from the fact that our reference triangle has area 1/2
 * and the determinant of the jacobian of an affine function mapping
 * an arbitrary triangle to the refernce triangle is |T|/|T_ref|
 */
template<class F, int Order>
double integrate(F&& f){
    double sum = 0.;
    for (int i = 0; i < Order; ++i) {
        auto [x,y] = Abscissas<Order>[i];
        sum += f(1 - x - y, x, y) * Weights<Order>[i];
    }
    return sum;
}

template<class F>
struct QuadratureRule {
    double (*integrator)(F&&);
    explicit QuadratureRule(int order){
        switch (order) {
            case 1:
            case 2:
            case 3 :
                integrator = integrate<F, 3>;
                break;
            case 4 :
                integrator = integrate<F, 4>;
                break;
            case 5 :
            case 6 :
                integrator = integrate<F, 6>;
                break;
            case 7 :
                integrator = integrate<F, 7>;
                break;
            case 8 :
            case 9 :
                integrator = integrate<F, 9>;
                break;
            case 10 :
            case 11 :
            case 12 :
                integrator = integrate<F, 12>;
                break;
            default :
                integrator = integrate<F, 13>;
                break;
        }
    }

    double operator()(F&& f){
        return integrator((F&&)f);
    }
};



