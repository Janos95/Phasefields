
#include "Viewer.h"
#include "ScopedTimer/ScopedTimer.h"

int main(int argc, char** argv){

    Phasefield::Viewer* viewer;
    {
        ScopedTimer t{"Application Start", true};
        viewer = new Phasefield::Viewer{{argc, argv}};
    }

    while(viewer->mainLoopIteration()) {
        if(viewer->isOptimizing) {
            viewer->runOptimization([&viewer]{ return viewer->mainLoopIteration(); });
            viewer->isOptimizing = false;
        }
    }
}

//MAGNUM_APPLICATION_MAIN(Phasefield::Viewer)