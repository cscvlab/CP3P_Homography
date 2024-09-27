# CP3P for Minimal Pose Estimation
This repository is the official implementation of the paper:

**A Concise P3P Method Based on Homography Decomposition**.

Concise Perspective-3-Point (CP3P) is a novel pose estimation method based on 2D homography decomposition. Under the minimal configuration consisted of three 3D-2D point correspondences, the homography between the object plane three 3D points lie on and the image is incomplete and has 2-DOF. Leveraging the proposed homography decomposition form, we obtain a fundenmental geometric relationship for the P3P problem, i.e., the similarity of two triangles, to construct quadratic constraints. This relationship can also be applied into early verification of ambiguous solutions with the fourth point correspondence, instead of the common post-processing reprojection. Moreover, the quadratic constraints and two unknown variables in most previous P3P methods are managed in a unified way. Last but no least, our derivation is extreamly concise and avoid unnecessary computations as much as possible (e.g., square root and division operations).  

**Links:** [[Project Page]](http://www.cscvlab.com/research/CP3P/)   

## Codes Explanation
### Project Structure Overview

| File/Directory name        |                                                                                                                                           Description                                                                                                                                            |  
|----------------------------|:------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------:|
| P3P_HD.h                   |                                                                            Our algorithm for solving the p3p problem, utilizing 4 and 3 point correspondences to compute a unique and multiple pose(s), respectively.                                                                             |
| demo_P3P.cpp               |                                                                                      A demo we provide for running all variants of our P3P methods, including 'Mul-4', 'Mul-3', 'Uni-3-1V', and 'Uni-3-3V'.                                                                                      |
| demo_Ransac.cpp            |                                             A demo we provided for running the RANSAC-based pose estimation process, allowing selection of different P3P solvers(Mul4, Mul3, Uni-3-1V, Uni-3-3V ), and optional use of EPNP and Ceres optimization.                                              |
| utils/functionSolver.h     | Methods for solving cubic equations based on a modified approach from Ding et. al's paper "Revisiting the P3P Problem" in CVPR2023 and quartic equations using techniques from Ke and Roumeliotis' paper "An efficient algebraic solution to the perspective-three-point problem" in CVPR2017. |           |
| utils/EPNP folder          |                             The folder contains functions for solving the Perspective-n-Point problem using the Efficient PnP (EPNP) algorithm from paper "EPnP: An Accurate O(n) Solution to the PnP Problem", which can be optionally used in the RANSAC process.                              |
| utils/colmap_util folder   |                                                      The folder contains utility functions from Colmap (https://github.com/colmap/colmap), used specifically for setting up and running Ceres optimization in the pose estimation process.                                                       |
| utils/cost_function folder |                                                                        The cost_function folder includes the objective functions for Ceres optimization, defining the error terms that guide the pose refinement process.                                                                        |
| Cambridge folder           |                 Several sets of data from the Cambridge KingsCollege scene, consisting of image points, corresponding 3D world points, and camera intrinsics, generated following the data preparation pipeline from the hloc(https://github.com/cvg/Hierarchical-Localization).                 |

## How to Run the Project
### 1. Install Dependencies
####   Install Eigen
Eigen is a C++ template library for linear algebra. You can install it from the source:

```bash
# Clone the Eigen source code
git clone https://gitlab.com/libeigen/eigen.git
cd eigen

# Create a build directory and install
mkdir build && cd build
cmake ..
sudo make install
```
#### Install Ceres Solver
Ceres Solver is a C++ library for solving non-linear least squares problems. You can install it with the following steps:

```bash
# Install Ceres dependencies
sudo apt-get install  liblapack-dev libsuitesparse-dev libcxsparse3 libgflags-dev libgoogle-glog-dev libgtest-dev

# Clone the Ceres source code
git clone https://ceres-solver.googlesource.com/ceres-solver
cd ceres-solver

# Create a build directory and install
mkdir build && cd build
cmake ..
make -j4
sudo make install
```
### 2. Configure and Build the Project
In the project’s root directory, run the following commands to generate the Makefile and build the project:

```bash
mkdir build
cd build
cmake ..
make
```
### 3. Run the Project
After building the project, you can run the executable with:

```bash
./build/CP3P_Homography
```

