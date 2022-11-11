## Linear Sparse Methods

Many imaging and compressed sensing applications seek sparse solutions to under-determined least-squares problems.
This directory contains sparse array optimization algorithms for direction of arrival estimation in ULA.
The following lines give particular information about each matlab file.

- [BPDN](BPDN.m): The from scratch implementation of Basis Pursuit Denoise algorithm for sparse array antenna recovery.
It is similar to Lasso implementation which is developed by the statistics community. However, the BPDN was developed by the
signal processing community. The Lasso and Basis Pursuit Denoising (BPDN) approaches of bounding the 1-norm of the solution have led to several computational algorithms.

- [Source_localization_V1](Source_localization_V1.m): It contains the MATLAb version implementation for matrix recovery. Then, the l1-norm objective function, iterated weighted l1-norm objective function, and rank reduction method are implemented via CVX.

- [SRACV](sracv.m): The implementation of the Sparse Reconstruction of Array Covariance Vector (SRACV) which is based on l1-norm optimization problem. Based on the sparse representation of array covariance vectors, a weighted L1-norm minimization is applied to the data model, in which the weighted vector can be obtained by taking advantage of the orthogonality between the
noise subspace and the signal subspace. By searching the sparsest coefficients of the array covariance vectors simultaneously, DOAs can be effectively estimated.

- [SRACV_MUSIC_CMP](sracv_music_cmp.m): It is the comparison of MUSIC and L1-SRACV algorithms with the RMSE (Root Mean Square Error) performance metric.

- [mySVD](mySVD.m): It is the accelerated singular value decomposition. Based on the size of X, mySVD computes the eigvectors of X*X^T or X^T*X first, and then convert them to the eigenvectors of the other.

- [SV](sv.m): It is based on the singular value decomposition (SVD). Then, its outputs are used for creating a temporary matrix in recovering signal.

- [L1-SVD](l1_svd.m): This is the l2,1 norm optimization for multiple snapshots for estimating dirction of arrival. It uses SV implementation in its own estimation.

- [Norm_lij](norm_lij.m): It is a specific version of norms to compute the norm of a matrix with i, j as its parameters.

- [RMSE](rmse.m): It is our used performance metric, which stands for Root Mean Square Error.

- [SLSA](SLSA.m): This implementation contaitns the following algorithms computation and comparison with the others:
    1) Single Sample Source locarization with BPDN in CVX
    2) l1-norm --- Multiple snapshots
    3) l2-1 norm optimization for multiple snapshots
    4) Beamformer - Capon - Music - l1-SVD

- [StoULA](StoULA.m): It contains some Gridless methods impelementation for estimating direction of arrival.
