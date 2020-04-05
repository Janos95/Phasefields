

#include "viewer.hpp"
#include "optimization_context.hpp"
#include "subdivision.hpp"
#include "brush.hpp"
#include "mesh_io.hpp"
#include "update_scene.hpp"
#include "primitives.hpp"


int main(int argc, char** argv) {

    Viewer viewer(argc, argv);
    Scene scene;
    auto* data = new PhasefieldData{}; //not sure how to pass this around...
    viewer.setScene(scene);
    viewer.insertEventCallbacks(
            UpdateScene{data},
            //MeshIO{scene},
            LoadPrimitives{data},
            //Brush{scene},
            Subdivision{data},
            //OptimizationContext{scene}
            );
    viewer.exec();
}
