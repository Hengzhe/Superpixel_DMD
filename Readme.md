## Multi-color superpixel method 
The code is for our paper published in Optics Express: [Multi-color complex spatial light modulation with a single digital micromirror device](https://opg.optica.org/oe/abstract.cfm?doi=10.1364/OE.494238)

We generalize the previous proposed superpixel method into a multi-color regime. Based on this generalized method, we demonstrated 3-color multi-plane holographic projection, and 2-color Airy beam generation. One can easily simulate those results with the code.

## Requirement
- MATLAB
- Toolbox: Statistics and Machine Learning Toolbox, Parallel Computing Toolbox

## Document
- "DLUT" includes the program to generate DLUT. I cannot upload my DLUT to github since they are too large. So, before you run those program, you need to re-generate the DLUT and KDTmodel. 
- "HolographicProjection" includes the code of simulation for holographic projection
- "MatDiffrac" include some code of Fresnel diffraction simulation.
- "util" includes some basic functions.
- "DMDforward.m" is a function for simulation of the light propagation from DMD plane to the output plane.
- "DMDgeneration_SP.m" is a function that generate DMDpattern from target light filed.
- "Demo.mlx" is a simple demonstration of our method, also saved as "Demo.pdf". 
You may follow the [demo](Demo.pdf) to learn how to use those code. 

## Citation
You may cite our paper: [Multi-color complex spatial light modulation with a single digital micromirror device](https://opg.optica.org/oe/abstract.cfm?doi=10.1364/OE.494238) if those code helps your work. 

Contact author: hz_yan@sjtu.edu.cn
