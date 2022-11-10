## Base Models

This file contains primary algorithms for direction of arrival estimation in uniform linear array antenna (ULA). The following lines give particular information about each file.

- [Covariance Matrix](Covariance.m): The implementation of covariance matrix for the output signal of ULAs.
The covariance matrix is the random process's second-order statistic, measured at the array sensors. It contains information about the sources in space (number, strength, direction) and can be used for source detection and separation. The number of independent spatial signals at the array's input and the R's rank is the same. It is worth saying algorithms based on covariance matrix inversion/decomposition suffer from source correlation. They can't be separated if the sources are highly correlated (have identical waveforms or too close a direction). And also, the system performance degrades in the highly correlated scenario.

- [Beamformere](Beamformer.m): The implementation of Beamformer method for the direction of arrival estimation task.
This estimator produces a single peak in the actual direction of arrival in a single target scene by measuring the power when directed in that direction. However, contributions from one source will bias the estimator output in other directions of arrival, even in an uncorrelated multi-source picture. As a result, the peaks move away from their actual arrival directions and toward one another. The peaks combine into a single peak, reducing resolution if both sources are situated inside the main beam.

- [Capon](Capon.m): The implementation of Capon method for the direction of arrival estimation task.
Beamforming is a spatial filtering technique used in signal processing. It is frequently utilized in sensor arrays to determine the signal's angle of arrival when noise is present. Capon Spectrum Beamforming is a minimum-variance, distortionless response technique. This means that it optimizes the beam's movement-to-noise ratio without affecting the signal's gain or phase.

- [MVDR](MVDR.m): The implementation of the Minimum Variance Distortionless Response (MVDR) for estimating direction of arrival.
The MVDR Beamformer is a data adaptive beamforming method whose objective is to reduce the output variance. The variance of the recorded signals is equal to the sum of the variances of the desired signal and the noise if the noise and the underlying desired signal are uncorrelated, which is often the case. The MVDR solution therefore aims to reduce this total in order to minimize the impact of the noise.

- [MUSIC](Music.m): The implementation of MUSIC(Multiple signal classification) algorithm for the direction of arrival estimation task. It is a subspace-based method for estimating the DOA of narrowband sources that arrive at a sensor array at the same frequency. When the number of components is known in advance, MUSIC outperforms simple approaches like identifying peaks of DFT spectra in the presence of noise because it uses this information to disregard the noise in its final report.

- [Root-MUSIC](Root_Music.m): The implementation of Root-MUSIC algorithm for the direction of arrival estimation task.
The root-MUSIC algorithm is a polynomial form of the MUSIC algorithm. This algorithm adopts the roots of a polynomial to replace the search for spatial spectrum in the MUSIC algorithm, reducing the calculation amount and improving estimation performance.

- [MUISIC-Like](Music_Like.m): The implementation of MUSIC-Like algorithm for the direction of arrival estimation task.
An algorithm called the MUSIC-like algorithm was initially proposed as an alternative method to the MUSIC algorithm for direction-of-arrival estimation. It was shown to have robust performance without requiring explicit model order estimation, particularly in low signal-to-noise ratio (SNR) scenarios.

- [ESPRIT](Esprit.m): The implementation of ESPRIT (Estimation of Signal Parameters via Rotational Invariance Techniques) algorithm for the direction of arrival estimation task.
The ESPRIT exploits the rotational invariance property of two identical subarrays and solves the eigenvalues of a matrix relating two signal subspaces. A simple way to construct identical subarrays is to select the first (M − 1) elements and the second to the Mth elements of a ULA.

- [Linear](linear.m): The implementation of Linear prediction model algorithm for the direction of arrival estimation task.
The linear prediction-based estimation procedure is commonly used in time series analysis for all pole modeling data. It has also been successfully used in array processing. In this case, one of the sensor outputs is predicted as a linear combination of the remaining sensor outputs at any instant. The predictor coefficients are selected to minimize the mean square error.

- [Demo](./Demo/): Contains demo files about implementations of the Base Models' algorithms.
-- Angle.m: Comparison of beamforming and capon method.
-- beamformer_capon_.m: Another version of beamforming and capon implementation
-- corr_beam_capon.m: Examining the beamforming and capon doa estimation methods with different implementation
-- flower.m: Implementation of estimating covariance matrix
-- music.m: Demo implementation for MUSIC algorithm
-- rmse_beam-capon.m: Computing RMSE for comparing capon and beamforming methods
-- shape.m: Estimating the direction of arrival capon and beamforming
