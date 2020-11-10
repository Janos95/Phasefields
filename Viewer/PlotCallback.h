//
// Created by janos on 10/6/20.
//

#pragma once

#include "Solver.h"
#include "Tree.h"
#include "RecursiveProblem.h"

#include <Corrade/Utility/FormatStl.h>

namespace Phasefield {

namespace Cr = Corrade;

//struct PlotCallback {
//
//    explicit PlotCallback(Solver::RecursiveProblem& pb, Node n) : problem{pb}, node(n), costs{pb.objectives.size()} {}
//
//    Solver::RecursiveProblem& problem;
//    Node node;
//    Array <Array<double>> costs;
//
//    Solver::Status::Value operator()(Solver::IterationSummary const&) {
//        for(auto& [f, hist, draw] : problem.objectives) {
//            double cost = 0;
//            f(node.phasefield(), node.temporary(), cost, nullptr, nullptr);
//            arrayAppend(hist, {float(hist.size()), float(cost)});
//        }
//        return Solver::Status::CONTINUE;
//    }
//
//    void generatePlot() {
//        std::string plot = "\\begin{tikzpicture}\n"
//                           "\\begin{axis}[\n"
//                           "height=9cm,\n"
//                           "width=9cm,\n"
//                           "grid=major,\n"
//                           "]\n";
//
//        size_t sampleCouunt = costs.front().size();
//        size_t step = Math::max(sampleCouunt/50, 1ul);
//        for(size_t k = 0; k < costs.size(); ++k) {
//            plot += "\n\\addplot coordinates {\n";
//            for(size_t j = 0; j < costs[k].size(); j += step) {
//                plot += Cr::Utility::formatString("({},{})\n", j, costs[k][j]);
//            }
//            plot += "};\n";
//            plot += Cr::Utility::formatString("\\addlegendentry{{ {} }}", f.functionalType));
//        }
//        plot += "\n\\end{axis}\n"
//                "\\end{tikzpicture}";
//
//        FILE *fp = fopen("/tmp/plot.tex", "w");
//        if (fp != nullptr) {
//            fputs(plot.data(), fp);
//            fclose(fp);
//        }
//    }
//};

}