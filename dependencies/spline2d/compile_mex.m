% Compile MEX files for the tspline package
mex -largeArrayDims spline2d_grad.c
mex -largeArrayDims spline2d_green.c