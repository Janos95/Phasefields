//
// Created by janos on 03.03.20.
//

#include "stopping_criteria.hpp"

#include <Corrade/TestSuite/Tester.h>
#include <random>
#include <Corrade/TestSuite/Compare/Numeric.h>

using namespace Corrade;

namespace {


    struct FW : TestSuite::Tester {

        FW()
        {
            addTests({&FW::positive, &FW::negative});
        }

        void positive(){
            std::random_device d;
            std::default_random_engine engine(d());
            std::uniform_real_distribution<double> dist(-1,1);
            for (int i = 0; i < 100; ++i) {
                auto x = dist(engine);
                if(0.85 < x && x < 0.95){
                    CORRADE_VERIFY(w1(x) > 0);
                    CORRADE_COMPARE(f1(x), 0);
                } else{
                    CORRADE_COMPARE(w1(x), 0);
                    CORRADE_VERIFY(f1(x)> 0);
                }
                auto eps = 1e-6;
                auto numericGradW = (w1(x+eps)-w1(x-eps))/(2*eps);
                CORRADE_COMPARE_WITH(numericGradW, w1g(x), TestSuite::Compare::Around{1e-6});
                auto numericGradF = (f1(x+eps)-f1(x-eps))/(2*eps);
                CORRADE_COMPARE_WITH(numericGradF, f1g(x), TestSuite::Compare::Around{1e-6});
            }
        }

        void negative(){
            std::random_device d;
            std::default_random_engine engine(d());
            std::uniform_real_distribution<double> dist(-1,1);
            for (int i = 0; i < 100; ++i) {
                auto x = dist(engine);
                if(-0.95 < x && x < -0.85){
                    CORRADE_VERIFY(w2(x) > 0);
                    CORRADE_COMPARE_WITH(f2(x), 0, TestSuite::Compare::Around{1e-6});
                } else{
                    CORRADE_COMPARE_WITH(w2(x), 0, TestSuite::Compare::Around{1e-6});
                    CORRADE_VERIFY(f2(x) > 0);
                }
                auto eps = 1e-6;
                auto numericGradW = (w2(x+eps)-w2(x-eps))/(2*eps);
                CORRADE_COMPARE_WITH(numericGradW, w2g(x), TestSuite::Compare::Around{1e-6});
                auto numericGradF = (f2(x+eps)-f2(x-eps))/(2*eps);
                CORRADE_COMPARE_WITH(numericGradF, f2g(x), TestSuite::Compare::Around{1e-6});
            }
        }

        W w1 = W(0.85,0.95);
        W w2 = W(-0.95,-0.85);
        F f1 = F(0.85,0.95);
        F f2 = F(-0.95,-0.85);

        WGrad w1g = WGrad(0.85,0.95);
        WGrad w2g = WGrad(-0.95,-0.85);
        FGrad f1g = FGrad(0.85,0.95);
        FGrad f2g = FGrad(-0.95,-0.85);

    };

}

CORRADE_TEST_MAIN(FW)