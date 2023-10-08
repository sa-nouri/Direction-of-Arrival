## Base Models for Direction of Arrival Estimation

This section contains fundamental algorithms for direction of arrival (DOA) estimation using uniform linear array antennas (ULAs). Each algorithm is briefly described below:

### [Covariance Matrix](Covariance.m)
- **Description**: Implements the covariance matrix for the output signal of ULAs.
- **Details**: The covariance matrix represents the second-order statistics of the random process, measured at the array sensors. It contains vital information about the sources in space, including their number, strength, and direction. Covariance matrices are often used for source detection and separation. However, algorithms based on covariance matrix inversion/decomposition can be limited by source correlation. They struggle to separate highly correlated sources or sources with nearly identical waveforms.

### [Beamformer](Beamformer.m)
- **Description**: Implements the Beamformer method for DOA estimation.
- **Details**: The Beamformer produces a single peak in the direction of the actual arrival of a single target source. However, in the presence of multiple sources, contributions from one source can bias the estimator's output in other directions, even if the sources are uncorrelated. This can lead to a reduction in resolution, with peaks moving away from their true arrival directions and potentially merging.

### [Capon](Capon.m)
- **Description**: Implements the Capon method for DOA estimation.
- **Details**: Capon Spectrum Beamforming is a spatial filtering technique used in signal processing, particularly in sensor arrays. It is employed to determine the angle of arrival of a signal when noise is present. The Capon method aims to optimize the beam's signal-to-noise ratio while preserving the signal's gain and phase.

### [MVDR](MVDR.m)
- **Description**: Implements the Minimum Variance Distortionless Response (MVDR) method for DOA estimation.
- **Details**: The MVDR Beamformer is a data-adaptive beamforming approach that seeks to minimize the output variance. It achieves this by reducing the impact of both the desired signal and noise. This method is effective when the noise and desired signal are uncorrelated.

### [MUSIC](Music.m)
- **Description**: Implements the MUSIC (Multiple Signal Classification) algorithm for DOA estimation.
- **Details**: MUSIC is a subspace-based method designed for estimating the DOA of narrowband sources arriving at a sensor array operating at the same frequency. It outperforms simple peak detection in the presence of noise, especially when the number of source components is known.

### [Root-MUSIC](Root_Music.m)
- **Description**: Implements the Root-MUSIC algorithm for DOA estimation.
- **Details**: Root-MUSIC is a polynomial form of the MUSIC algorithm. It employs the roots of a polynomial to identify spatial spectrum peaks, reducing computational complexity while improving estimation performance.

### [MUSIC-Like](Music_Like.m)
- **Description**: Implements the MUSIC-Like algorithm for DOA estimation.
- **Details**: The MUSIC-Like algorithm serves as an alternative to the MUSIC algorithm for DOA estimation. It exhibits robust performance, particularly in low signal-to-noise ratio (SNR) scenarios, without requiring explicit model order estimation.

### [ESPRIT](Esprit.m)
- **Description**: Implements the ESPRIT (Estimation of Signal Parameters via Rotational Invariance Techniques) algorithm for DOA estimation.
- **Details**: ESPRIT leverages the rotational invariance property of two identical subarrays to estimate eigenvalues of a matrix related to signal subspaces. Identical subarrays are constructed by selecting specific elements from a ULA.

### [Linear Prediction](linear.m)
- **Description**: Implements a linear prediction model algorithm for DOA estimation.
- **Details**: Linear prediction-based estimation is commonly used in time series analysis and has been successfully applied in array processing. This method predicts one sensor's output as a linear combination of the remaining sensor outputs, minimizing mean square error.

### [Demo](./Demo/)
This directory contains demo files related to the implementation of the Base Models' algorithms. These files include:

- **Angle.m**: A comparison of beamforming and Capon methods.
- **beamformer_capon_.m**: An alternative implementation of beamforming and Capon methods.
- **corr_beam_capon.m**: Examination of beamforming and Capon DOA estimation methods with different implementations.
- **flower.m**: Implementation of covariance matrix estimation.
- **music.m**: A demonstration of the MUSIC algorithm.
- **rmse_beam-capon.m**: Computation of the root mean square error (RMSE) for comparing Capon and beamforming methods.
- **shape.m**: Estimation of DOA using Capon and beamforming.

Please explore these resources to gain insights into DOA estimation techniques and array signal processing. If you have any questions or require further assistance, feel free to reach out for clarification or discussion.
