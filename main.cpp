
#include "Viewer.h"

int main(int argc, char** argv){
    Phasefield::Viewer viewer{argc, argv};

    while(viewer.mainLoopIteration()){
        if(viewer.isOptimizing){
            auto terminationType = viewer.runOptimization([&viewer]{ return viewer.mainLoopIteration(); });
            if(terminationType == Phasefield::Solver::Status::ABORT)
                Magnum::Debug{} << "User terminated optimization";
            viewer.isOptimizing = false;
        }
    }
}