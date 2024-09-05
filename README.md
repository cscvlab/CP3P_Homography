# CP3P for Minimal Pose Estimation
This repository is the official implementation of the paper:

**A Concise P3P Method Based on Homography Decomposition**.

Concise Perspective-3-Point (CP3P) is a novel pose estimation method based on 2D homography decomposition. Under the minimal configuration consisted of three 3D-2D point correspondences, the homography between the object plane three 3D points lie on and the image is incomplete and has 2-DOF. Leveraging the proposed homography decomposition form, we obtain a fundenmental geometric relationship for the P3P problem, i.e., the similarity of two triangles. This property can also be applied into early verification with the fourth point correspondence, instead of the common post-processing reprojection. Moreover, the two unknown variables in most previous P3P methods are managed in a unified way. Last but no least, our derivation is extreamly concise and avoid unnecessary computations as much as possible (e.g., square root and division operations).  

**Links:** [[Project Page]](http://www.cscvlab.com/research/CP3P/)   

## Main Results
FLOPs

Accuracy w.o. noise

Time

RANSAC

## Codes Explanation
Files

Run


