//
// Created by janos on 03.03.20.
//

#include "dijkstra.hpp"
#include "bfs.hpp"

#include <Corrade/TestSuite/Tester.h>
#include <random>
using namespace Corrade;

namespace {


    struct DijkstraTest : TestSuite::Tester {

        using G = std::vector<std::vector<std::pair<int, double>>>;
        G g1 = {
                {{1, 5}, {2, 1}},
                {{3, 6}, {4, 9}},
                {{3, 10}, {4, 2}},
                {{5, 8}, {6, 7}},
                {{5, 3}},
                {{6, 4}},
                {},
        };


        G g2 = {
            {{1,5}, {2,1}},
            {{3, 6}, {4,9}, {0, 5}},
            {{3, 10}, {4,2}, {0, 1}},
            {{5, 8}, {6,7}, {1, 6}, {2, 10}},
            {{5, 3}, {1, 9}, {2, 2}},
            {{6, 4}, {3, 8}, {4, 3}},
            {{3, 7}, {5, 4}},
        };

        G g3 = {
                {{1,1}, {2,1}},
                {{3, 1}, {4,1}},
                {{3, 1}, {4,1}},
                {{5, 1}, {6,1}},
                {{5, 1}},
                {{6, 1}},
                {},
        };


        DijkstraTest()
        {
            addTests({
                &DijkstraTest::compareToBfs,
                &DijkstraTest::weighted,
                &DijkstraTest::weightedFull,
                &DijkstraTest::earlyStopping});
        }

        void weighted(){
            Dijkstra dijk(g1);
            dijk.run(0, {});
            auto r = dijk.getShortestPathReversed(0, 6);
            auto l = std::accumulate(r.begin(), r.end(), 0.
                    ,[&](const auto s, const graph::Edge& e)
                    {
                        auto it = std::find_if(g1[e.a].begin(), g1[e.a].end(), [&](auto& p){ return p.first == e.b; });
                        return s + it->second;
                    });
            CORRADE_COMPARE(l, 10.);
        }


        void weightedFull(){
            Dijkstra dijk(g2);
            dijk.run(0, {});
            auto r = dijk.getShortestPathReversed(0, 6);
            auto l = std::accumulate(r.begin(), r.end(), 0.
                    ,[&](const auto s, const auto& e)
                     {
                         auto it = std::find_if(g2[e.b].begin(), g2[e.b].end(), [&](auto& p){ return p.first == e.a; });
                         return s + it->second;
                     });
            CORRADE_COMPARE(l, 10.);
        }


        void earlyStopping(){
            Dijkstra dijk(g2);
            auto stop = [](auto n){ return n == 4; };
            dijk.run(0, {stop});
            auto r = dijk.getShortestPathReversed(0, 4);
            auto l = std::accumulate(r.begin(), r.end(), 0.
                ,[&](const auto s, const auto e)
                 {
                     auto it = std::find_if(g2[e.b].begin(), g2[e.b].end(), [&](auto& p){ return p.first == e.a; });
                     return s + it->second;
                 });
            CORRADE_COMPARE(l, 3.);
        }

        void compareToBfs(){
            Dijkstra dijk(g3);
            BreadthFirstSearch bfs(g3);
            dijk.run(0, {});
            bfs.run(0);
            for (int i = 1; i < 7; ++i) {
                auto r1 = dijk.getShortestPathReversed(0, i);
                auto r2 = bfs.getShortestPathReversed(0, i);
                auto eq = std::equal(r1.begin(), r1.end(), r2.begin());
                CORRADE_VERIFY(eq);
            }
        }
    };

}

CORRADE_TEST_MAIN(DijkstraTest)