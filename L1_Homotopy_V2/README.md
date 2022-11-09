## L1 Homotopy Package

l1homotopy is a highly versatile homotopy program that can solve a variety of L1-norm minimization problems using a warm start.

l1homotopy.m is the main function that solves the following homotopy program:

minimize_x  \|W x\|_1 + 1/2*\|Ax-y\|_2^2 + (1-epsilon)u'x,

u is defined as u = -W*sign(xh_old)-A'*(A*xh_old-y)
xh_old is an arbitrary warm-start vector (or a zero vector if no warm-start is available)

Scripts for different problems are also included in this package to demonstrate the use of l1homotopy:

- demo_BPDN -- solves LASSO/BPDN problem from scratch
- demo_posBPDN -- solves BPDN problem with positive sign constraint on the estimate
- demo_dynamicX -- updates the solution for BPDN as the signal changes
- demo_dynamicSeq -- updates the signal as sequential measurements are added
- demo_rwtL1-- solves iterative reweighting for BPDN
- demo_dynamicRWT -- iteratively updates weights in the L1-norm cost while estimating a time-varying signal
- demo_streamingLOT -- iteratively estimates a streaming signal using lapped orthogonal transform as the representation basis
- demo_KalmanRWT -- iteratively estimates a streaming signal that follows a linear dynamic model
- and ...

You may need to compile mex codes for
1. matrix-vector product of the form A_Gamma x_Gamma and A_Gamma^T y
2. realnoiselet
3. Wavelets

See compile.m for further details.
Add all the folders in MATLAB path or only those that are required for each solver.

Other than L1 decoding and adaptive reweighting methods, these homotopy programs can also be solved using l1homotopy.m (for the LASSO/BPDN formulation).

### Pursuits_Homotopy (Standard homotopy solvers

- Basis pursuit denoising (BPDN) homotopy
- BPDN_homotopy_function.m
- BPDN_homotopy_demo.m

- Dantzig selector (DS) homotopy based on primal-dual pursuit
- DS_homotopy_function.m
- DS_homotopy_demo.m

### DynamicX (Homotopy update for time-varying sparse signals)

- BPDN
- DynamicX_BPDN_function.m
- DynamicX_BPDN_demo.m
- DynamicX_BPDN_Visual_demo.m

- DS
- DynamicX_DS_function.m
- DynamicX_DS_demo.m

- Simulations (use functions from GPSR and FPC_AS)
- Simulation_DynamicX_BPDN.m
- Simulation_DynamicX_BPDN_Pathological.m
- Simulation_DynamicX_BPDN_Wavelet.m % This script uses WaveLab functions.

### DynamicSeq (Homotopy update for sequential measurements)

- BPDN
- DynamicSeq_BPDN_function.m
- DynamicSeq_BPDN_demo.m

- DS
- DynamicSeq_DS_function.m
- DynamicSeq_DS_demo.m

- Simulations (use functions from GPSR and FPC_AS)
- Simulation_DynamicSeq_BPDN.m

### Decoding (Correct sparse errors in an encoded signals)

- L1 decoding
- l1Decode_homotopy_fast.m
- l1Decode_homotopy_qr.m
- Simulation_l1Decode.m

- Robust error correction
- DynamicSeq_REC_function.m
- DynamicSeq_REC_demo.m
- Simulation_DynamicSeq_REC.m

### Iterative and adaptive reweighting (WeightedBPDN folder)

More details are availabe in the [WeightedBPDN-README](./WeightedBPDN/README.md).