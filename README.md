This repository contains the numerical code for my masterthesis.
It provides a gui application for the hierarchical optimization of phasefields,
for computing surface segmentations which can be conformally mapped
to the plane with low area distortion (I know, its a mouth full).

There is also an implementation of a brush, with which one can edit the phasefield
and various visualization options.

**Introduction**

You can try it out in your [broiwser](https://janos95.gitlab.io/wasm-test/), but be aware the 
webasm build is pretty buggy atm.

Here is a teasure of the optimization result for varying interface length penalties.
![](images/image.png)

**building**

This repo uses cmake for build file generation. 
If you have Eigen installed on your system and didn't forget
to init the git submodules you should be able to build the 
project using e.g. ninja as follows:

```bash
mkdir build && cd build && cmake .. -GNinja && ninja
```

There are some more dependancies which you can enable using cmake options in brackets :
- Loading meshes using Assimp (PHASEFIELD_WITH_IO).
- Exporting videos via ffmpeg (PHASEFIELD_WITH_IO)
- Solving constrained optimization problems using  Ipopt (PHASEFIELD_WITH_IPOPT).
- Checking derivatives via AD using ADOL-C (PHASEFIELD_WITH_ADOLC).
- Solving the yamabe equation faster using Suite Sparse (PHASEFIELD_WITH_SUITESPARSE).
- A different lBFGS implementation provided by ceres-solver (PHASEFIELD_WITH_CERES).
- Threading using TBB (PHASEFIELD_WITH_TBB).


