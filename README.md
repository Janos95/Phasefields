This repository contains the numerical code for my masterthesis.

**Introduction**

You can try it out [youself](phasefields.org)!

The main target is a gui application which can be used to perform phasefield optimizations
wrt. different objectives. Teaser image 

![](images/image.png)

**building**

This repo uses cmake for build file generation. 
To build using ninja : 

```bash
mkdir build && cd build && cmake .. -GNinja && ninja
```

