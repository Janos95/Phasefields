

#include "viewer.hpp"
#include "optimization_context.hpp"
#include "subdivision.hpp"
#include "brush.hpp"
#include "mesh_io.hpp"
#include "primitives.hpp"
#include "colormaps.hpp"
#include "shader_options.hpp"

/*
 * @todo in release: load u, change shader -> sigsev
 */
int main(int argc, char** argv) {

    Viewer viewer(argc, argv);
    auto data = new PhasefieldData; //@todo leaks
    //ownership taken by viewer
    viewer.insertEventCallbacks({
            //MeshIO{scene},
            new LoadPrimitives{*data},
            new Brush{*data},
            new Subdivision{*data},
            new OptimizationContext{*data},
            new ColorMap{*data},
            new ShaderOptions{*data}
    });
    viewer.exec();
}
