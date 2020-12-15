***Introduction***

This repository contains the numerical code for my master thesis.
It provides a GUI application for the hierarchical optimization of phasefields.
Some of the highlights are :

- A phasefield parameterization of the Yamabe energy,
which can be used to compute surface segmentations which can be conformally mapped
to the plane with low area distortion.

- A topological penalty term which enforces connectedness of the segments.

- A hierarchical segmentation scheme, which combines multiple phasefields, to 
compute segmentations with more than two segments.

The GUI also support some other niceties like editing phasefields
using a brush-like tool, mesh IO via Assimp, making and exporting videos
using ffmpeg, color optimization for segmentations etc. 

You can try it out in your [browser](https://janos95.gitlab.io/wasm-test/), but be aware the 
webasm build is pretty buggy atm. Also its using blender style navigation, so dont forget to press
the middle mouse button when dragging.

Here is a teaser image of the optimization result for varying interface length penalties.
![](images/image.png)

***Building***

Before building make sure that git submodules are downloaded.
You can download them via:

```bash
git submodules update --init --recursive
```

Before compiling you also need to make sure you have Eigen 
and SuiteSparse installed. If not already present I recommend
installing them using you package manager.

Build files are generated using cmake. 
If you want to build the project using e.g.
ninja, you can use the following command:

```bash
mkdir build && cd build && cmake .. -GNinja -DCMAKE_BUILD_TYPE=Release && ninja
```

There are some more dependencies which you can enable using the cmake options in paranthesis :
- Loading meshes using Assimp (PHASEFIELD_WITH_IO).
- Exporting videos via ffmpeg (PHASEFIELD_WITH_IO)
- Solving constrained optimization problems using  Ipopt (PHASEFIELD_WITH_IPOPT).
- Checking derivatives via AD using ADOL-C (PHASEFIELD_WITH_ADOLC).
- Solving the yamabe equation faster using Suite Sparse (PHASEFIELD_WITH_SUITESPARSE).
- A different lBFGS implementation provided by ceres-solver (PHASEFIELD_WITH_CERES).
- Threading using TBB (PHASEFIELD_WITH_TBB).

***Related Work***

The diffuse Yamabe energy draws inspiration from the paper 
'Variational Surface Cutting' for which there is also code 
available [here](https://github.com/nmwsharp/variational-surface-cutting)


