# LCD-Module-Process-Mura-Defect-Detection
 Use Laplacian Of Gaussian function filter to generate 3 types of convolution and detect 3 types of mura respectively. The first type: 2D 3x3 convolution to detect defect point and line mura. Type 2: 1D 1x15 convolution, detect H-Band mura The third type: 1D 15x1 convolution, detect V-Band mura
# Laplacian Of Gaussian(LOG) Filter principle
The principle of detecting defects is mainly based on the LOG filter for image processing. The LOG filter is mainly combined with the second-order Laplacian
A filter with differential zero crossing point detection edge function and Gaussian noise smoothing function
