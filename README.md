# Light Field Inpainting Propagation via Low Rank Matrix Completion

## Description
This code implements the method described in the paper:

Mikaël Le Pendu, Xiaoran Jiang and Christine Guillemot, ["Light Field Inpainting Propagation via Low Rank Matrix Completion"](https://hal.archives-ouvertes.fr/hal-01720622/file/LFInpaintTIP.pdf), IEEE Transactions on Image Processing, vol. 27, no. 4, pp. 1981–1993, Apr. 2018.

Please cite this paper when using this code in your research.
More details and supplementary materials can also be found on [the project's webpage](https://www.irisa.fr/temics/demos/lightField/InpaintMC/LFinpaintMC.html)

## Content
- LF_utilities

Basic tools for light fields(display, refocus, save as a png sequence).

- homography_tools

Tools for generating random homographies and warping.

- LRTC

Low rank matrix completion code based on ADMM method. Both trace norm and rank minimization are available.

- misc_utilities

Additional tools.

- LFData

Directory for input data. A Light Field example is given.

- readHCI_LF.m

Load light field data from png files as well as mask and segmentation data.

- RandomHomography.m

Main inpainting propagation program.

## Usage
1. Launch readHCI_LF to read a light field with pre-inpainted central view + Mask.

(Original light field stored in M_org matrix)

2. Launch RandomHomography for propagation to all views.

(Inpainted light field stored in X_rec matrix)

3. Display (or save) lightField matrix with displayLFMatrix (or saveLFMatrix)

e.g. Video = displayLFMatrix(X_rec,imgSize,nU,nV,1);
