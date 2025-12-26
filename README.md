# 524FinalProject: Dynamic Deformables

https://github.com/user-attachments/assets/0b72c7e1-f630-4530-8afb-82ff3ead77ce

*NOTE: This is a work-in-progress that I am hoping to add to in the future!*

This repository contains my implementation of the isotropic volumetric deformation portion of Dynamic Deformables by Theodore Kim and David Eberle (2022) as my final project for CPSC524: Geometric Modeling. It specifically implements the BDF-1 integration scheme for tetrahedral meshes, with support for Rayleigh Damping and external forces. The elastic energies implemented thus far are As-Rigid-As-Possible (ARAP) and Stable Neo-Hookean (SNH). Because there is a reference implementation (HOBAK), I admit that took heavy inspiration from HOBAK's organizational structure. However, almost all of the code is self written. See the **Academic Honesty** section below for more details on this. Links to the references used are below:

- Dynamic Deformables: https://www.tkim.graphics/DYNAMIC_DEFORMABLES/
- HOBAK: https://github.com/theodorekim/HOBAKv1

# Implementation Notes
This project was built primarily using Eigen [Guennebaud, 2013] for matrix computation and Polyscope [Sharp, 2019] for interactive visualization.

- Eigen: https://gitlab.com/libeigen/eigen
- Polyscope: https://github.com/nmwsharp/polyscope

# Building the Project
*NOTE: This project has only been tested on Apple Silicon (M1). For other platforms, you may need to modify the CMakelists.txt*

Begin by cloning the repository:
```
git clone --recurse-submodules git@github.com:dcui-03/524FinalProject.git
cd ./524FinalProject
```
Create a deps folder and clone Polyscope and libigl:
```
mkdir deps && cd ./deps
git clone --recurse-submodules https://github.com/libigl/libigl.git
git clone --recurse-submodules https://github.com/nmwsharp/polyscope.git
```
Create a build folder and build the project:
```
cd .. && mkdir build && cd ./build
cmake ..
make -j4
```
This will create an executable in the build folder called *deformation_app*.

# Generating a Deformation

To run the project, you will need a path to a TOBJ file (see Data folder or HOBAK (https://github.com/theodorekim/HOBAKv1)) as well as a *generic path* to your output location. Since this project performs offline deformation on potentially large meshes, I opted to save deformations as sequences of files which are loaded as needed on a per-frame basis during playback. This requires a generic output path (ex. ./output/mymesh.txt) which describes the folder and file prefix to dump sequences to. The general format to run the project from the build folder is as follows:

```
./deformation_app <path_to_input_tobj> <generic_output_path>
```

This will launch a Polyscope window with your input tet mesh. There are several GUI buttons which allow you to control simulation parameters, such as the integration scheme (BDF-1 or BDF-2), the Elastic Energy being used (SNH or ARAP), as well as the simulation timespan (in frames and frames per second) and material parameters (Young's Modulus and Poisson's Ratio). You may also use this to toggle Rayleigh Damping and an External Force direction.

Another button is called *Apply Constraint*. Click the button to activate Constraint Mode, then select the vertices you wish to constrain. After that, click the button again to exit Constraint Mode, or press *Reset Constraints* to clear all of the currently selected constraints. Currently, applying a constraint will fully constrain the position of the selected vertex in all directions.

When you are ready to simulate, press *Compute Deformation*. This will compute the deformation sequence for the specified parameters. The simulator will dump the output deformation sequence to the specified output path.

At this point, you may press *Playback* to visualize the animation. If desired, you can pause the playback to change the frame rate using the *Steps Per Second* slider. You may instead choose to toggle the slider by pressing *Scrubbing*, which allows you to load per frame data as desired.

# Deformation Playback

If you have already generated a deformation sequence but have closed the Polyscope window, you have the option to play it back using the dumped files. To do so, you will need to run Playback Mode:

```
./deformation_app <path_to_input_tobj> <generic_output_path> --p <index_of_first_frame> <number_of_frames_to_visualize>
```

For example, if I wanted to visualize 75 frames of a saved simulation sequence beginning from frame 0, I might run the following:
```
./deformation_app ../data/bunny_6.tobj ../output/bunny_6.txt --p 0 75
```

This will again open the Polyscope window, but you will be locked to Playback Mode, which only allows you to play back or scrub the specified sequence.

# Future TODOs
- Anisotropic energies
- Collisions and contact forces
- 2D triangle meshes
- 3D thin shells
- Load tet meshes using libigl::tetgen
- Direction-independent application of constraints
- Find a better file saving/loading format...

# Academic Honesty
Any functions taken or heavily inspired by an external (non-AI) source is explicitly labeled as such with comments. These are primarily from HOBAK itself, with the functions in question being in /src/utils/MatrixUtils:
- computedFdx(): Copied from Dynamic Deformables [Kim & Eberle, 2022] (MATLAB code)
- partialJpartialF(): Directly copied from HOBAK
- buildTwistAndFlipEigenvectors(): Heavily inspired by HOBAK
- buildScalingEigenvectors(): Heavily inspired by HOBAK

Additionally, all functions for reading and writing files in /src/utils/IOutils are written largely by AI. All other code was written by me :)
