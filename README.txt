SANDBOXPIV SOFTWARE SUITE v0.11 --------------------------------------------------------------------

This PIV implementation was designed specifically for use with the Yalebox analog model of crustal deformation.  Features are:

- Direct spatial doamin computation of the correlation function, allowing for large displacements of small sampling windows.

- Iterative grid refinement, allowing for accurate results at high resolution

- Flexible vector validation

- Transformation from pixel to world coordinates

- Smoothing and gradient functions that ignore values outside the sand, eliminating spurious edge gradients

Enjoy!

WORKFLOW --------------------------------------------------------------------------------------------

1) Begin with a suitable image series (corrected lens distortion and perspective, masked regions outside the sand with 0's, rotated such that the tabletop is horizontal).
2) Explore settings either manually or using testparams.m to determine a suitable combination for the particular image series at hand.
3) Define the run paramenters using writeinput.m
4) Call runpivseries.m to run sandboxpiv.m for each image pair in the series 
5) Define the world coordinate system with getpivwoco.m and a suitable world coordinate image
6) Convert all piv data files to world coordinates with datapostprocess.m
7) Manipulate as you please, using nangrad2.m and nansmooth2.m as needed to reduce spurious edge gradients.
