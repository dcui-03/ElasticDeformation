# Nonlinear Elastic Deformation With Dynamic Deformables

https://github.com/user-attachments/assets/0b72c7e1-f630-4530-8afb-82ff3ead77ce

This repository contains my implementation and extensions of the isotropic volumetric deformation portion of Dynamic Deformables by Theodore Kim and David Eberle (2022) as my final project for UBC CPSC524: Geometric Modeling. It mainly implements the BDF-1 and BDF-2 integration scheme for tetrahedral meshes, with support for Rayleigh Damping and external forces. The elastic energies implemented thus far are As-Rigid-As-Possible (ARAP) and Stable Neo-Hookean (SNH). Because there is a reference implementation (HOBAK), I admit that took heavy inspiration from HOBAK's organizational structure. However, almost all of the code is self written. See the **Academic Honesty** section below for more details on this. Links to the references used are below:

- Dynamic Deformables: https://www.tkim.graphics/DYNAMIC_DEFORMABLES/
- HOBAK: https://github.com/theodorekim/HOBAKv1

Moving beyond Dynamic Deformables' foundation, I have also added support for quasi-static deformation (code needs some cleaning) via Newton iterations with line search. Note this is not very robust yet as I'm still fiddling with it. Lastly, I added support for abs-style projection when the Poisson's Ratio is close to 0.5 to stabilize the system, as recommended by the following paper:

- *Stabler Neo-Hookean Simulation: Absolute Eigenvalue Filtering for Projected Newton* by Chen et al. (2022) (https://www.cs.columbia.edu/cg/abs-psd/)

# Implementation Notes
This project was built primarily using Eigen [Guennebaud, 2013] for matrix computation and Polyscope [Sharp, 2019] for interactive visualization.

- Eigen: https://gitlab.com/libeigen/eigen
- Polyscope: https://github.com/nmwsharp/polyscope

# Building the Project
*NOTE: This project has only been tested on Apple Silicon (M1). For other platforms, you may need to modify the CMakelists.txt*

Begin by cloning the repository:
```
git clone --recurse-submodules git@github.com:dcui-03/ElasticDeformation.git
cd ./ElasticDeformation
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

To run the project, you will need a path to a TOBJ file (see Data folder or HOBAK (https://github.com/theodorekim/HOBAKv1)) as well as a *generic path* to your output location. Since this project performs offline deformation on potentially large meshes, I opted to save deformations as sequences of files which are loaded as needed on a per-frame basis during playback. This requires a generic output path (ex. ./output/mymesh.txt) as a TXT which describes the folder and file prefix to dump sequences to. The general format to run the project from the build folder is as follows:

```
./deformation_app <path_to_input_tobj> <generic_output_path>
```

This will launch a Polyscope window with your input tet mesh. There are several GUI buttons which allow you to control simulation parameters, such as the integration scheme (Quasi-static, BDF-1, or BDF-2), the Elastic Energy being used (SNH or ARAP), as well as the simulation timespan (in frames and frames per second) and material parameters (Young's Modulus and Poisson's Ratio). You may also use this to toggle Rayleigh Damping and an External Force direction.

You can also apply constraints at vertices using the *Apply Constraints* option in conjunction with the specified axes, as well as move vertices around using a Polyscope gizmo by pressing *Move Vertex*.

When you are ready to simulate, press *Compute Deformation*. If doing time integration, this will compute the deformation sequence for the specified parameters. The quasistatic deformation treats the result as a single frame (may change this in the future). The simulator will dump the output deformation sequence to the specified output path.

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
- Static Deformer (ex. via Laplacian) for posing

# Academic Honesty
Any functions taken or heavily inspired by an external (non-AI) source is explicitly labeled as such with comments. These are primarily from HOBAK itself and are labeled as such in the code. Much thanks to the authors of Dynamic Deformables for their clear explanations and helpful code bits.

Additionally, all functions for reading and writing files in /src/utils/IOutils are written largely by AI. All other code was written by me :)
