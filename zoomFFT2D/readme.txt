1D propagation:
---------------
Fraunhofer and Fresnel propagation routines for radially symmetrical fields, and a few supporting routines.

Run "tests.m" for an example of usage.

zoomFFT2D (Fraunhofer propagation):
-----------------------------------
Includes 4 versions of the zoomFFT routines, which is a fast implementation of a 2D Fourier Transform on arbitrarily specified input and output periodic grids. See the header of each function for detailed descriptions.

zoomFFT_realunits.m is the Fraunhofer integral (to within a constant quadratic phase factor).

czt.m is a subroutine on which the zoomFFT routines are based. It is a standard function included in the MATLAB signal processing toolkit.

czt_memory.m is a (non-standard) version of czt that is slower but more memory efficient.

Fresnel2D:
---------------
Two implementations of Fresnel propagation, appropriate in different circumstances. See "FresnelPropagateTest.m" for details on usage.