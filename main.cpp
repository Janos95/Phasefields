

#include "viewer.hpp"
#include "optimization_context.hpp"

int main(int argc, char** argv) {

    Viewer viewer(argc, argv);

    OptimizationContext optimizationContext;
    viewer.tickCallbacks.emplace_back([&](Scene*& s){ optimizationContext.updateScene(s); });
    viewer.menuCallbacks.emplace_back([&](auto& c){ optimizationContext.showMenu(c); });
    viewer.exec();
}
