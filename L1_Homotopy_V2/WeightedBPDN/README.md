## Weighted BPDN

This package contains MATLAB code for solving iterative reweighted and adaptive reweighted L1 problems using the homotopy methods described in the paper:

"Fast and accurate algorithms for re-weighted L1 norm minimization," by M. Salman Asif and Justin Romberg School of ECE, Georgia Tech.

For usage details, consult demo_wtBPDN.m and demo_adpWBPDN.m

This package also includes the following algorithms by
other authors (to allow running the comparative tests):

The SpaRSA algorithm, which can be downloaded from
http://www.lx.it.pt/~mtf/SpaRSA/

The YALL1 algorithm, which can be downloaded from
http://yall1.blogs.rice.edu/

The SPGL1 algorithm, which can be downloaded from
http://www.cs.ubc.ca/~mpf/spgl1/

This code is in development stage; any comments or bug reports are very welcome.

To reproduce results in the paper, use the following scripts:

job_wtBPDN_WAVE.m (Blocks and HeaviSine)

job_adpWBPDN.m (grayscale images)

May need to compile some mex files for matrix-vector product noiselet and wavelet transforms:  See compile.m in the main directory.

The package also includes a modified version of SpaRSA (renamed as SpaRSA_adpW)
that adaptively changes weights at every continuation step.
adp_wt = 1 sets the flag for adaptive reweighting

See SpaRSA_adpW.m for the details.
