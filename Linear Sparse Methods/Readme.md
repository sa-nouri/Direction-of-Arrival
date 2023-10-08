## Linear Sparse Methods

In numerous imaging and compressed sensing applications, the pursuit of sparse solutions to under-determined least-squares problems is paramount. This directory houses a collection of sparse array optimization algorithms designed for direction of arrival estimation in Uniform Linear Array (ULA) configurations. Below, you'll find detailed information about each MATLAB file:

- [BPDN (BPDN.m)]: This file presents a ground-up implementation of the Basis Pursuit Denoise (BPDN) algorithm for recovering sparse array antenna solutions. BPDN is akin to the Lasso implementation developed by the statistics community. However, it was specifically developed by the signal processing community. Both Lasso and Basis Pursuit Denoising (BPDN) approaches, which bound the 1-norm of the solution, have paved the way for several computational algorithms.

- [Source Localization Version 1 (Source_localization_V1.m)]: This file contains MATLAB implementations for matrix recovery. It includes the l1-norm objective function, iterated weighted l1-norm objective function, and rank reduction method, all implemented via CVX.

- [SRACV (sracv.m)]: Here, you'll find the implementation of the Sparse Reconstruction of Array Covariance Vector (SRACV) algorithm, based on l1-norm optimization principles. Leveraging the sparse representation of array covariance vectors, this approach employs weighted L1-norm minimization within the data model. The weighted vector is obtained by capitalizing on the orthogonality between the noise subspace and the signal subspace. By simultaneously identifying the sparsest coefficients of the array covariance vectors, it effectively estimates directions of arrival (DOAs).

- [SRACV_MUSIC_CMP (sracv_music_cmp.m)]: This file offers a comparison between the MUSIC and L1-SRACV algorithms using the Root Mean Square Error (RMSE) performance metric.

- [mySVD (mySVD.m)]: This accelerated singular value decomposition (SVD) method adapts its computations based on the size of the input matrix X. It computes the eigenvectors of either X*X^T or X^T*X first and then converts them into eigenvectors of the other matrix.

- [SV (sv.m)]: This implementation is based on singular value decomposition (SVD). Its outputs contribute to the creation of a temporary matrix used in signal recovery.

- [L1-SVD (l1_svd.m)]: This algorithm optimizes the l2,1 norm for multiple snapshots to estimate the direction of arrival. It utilizes the SV implementation within its estimation process.

- [Norm_lij (norm_lij.m)]: This file provides specific versions of norm calculations to compute the norm of a matrix with parameters i and j.

- [RMSE (rmse.m)]: RMSE stands for Root Mean Square Error and serves as the performance metric used in our evaluations.

- [SLSA (SLSA.m)]: This implementation encompasses the computation and comparison of various algorithms, including Single Sample Source Localization with BPDN in CVX, l1-norm for multiple snapshots, l2-1 norm optimization for multiple snapshots, Beamformer, Capon, Music, and l1-SVD.

- [StoULA (StoULA.m)]: This file contains implementations of Gridless methods for estimating the direction of arrival in sparse configurations.

These resources provide an in-depth exploration of linear sparse methods for direction of arrival estimation in ULAs. If you have any questions, require further clarification, or wish to discuss any aspects of these methods, please feel free to reach out.
