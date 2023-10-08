## Circular Array Antenna

The Uniform Circular Array (UCA) features elements uniformly distributed along a circle with a radius, R, where the origin of the coordinate system is placed at the array's center. Similar to other arrays, the UCA may also be susceptible to grating lobes. To mitigate this, it is common practice to select the radius in such a way that the minimum element separation closely aligns with λ/2, similar to the case of a Uniform Linear Array (ULA).

Here is specific information about each implemented MATLAB file for circular array antennas:

- [Circular-MUSIC (CircMusicDavis.m)]: This is a two-dimensional Multiple Signal Classification (MUSIC) method tailored for UCA.

- [Circular-Vandermond (Circular_vandermond.m)]: This approach employs Vandermonde decomposition for estimating the direction of arrival in UCA.

- [Circular-MUSIC-Davis (CircMusicDavis.m)]: This version of the MUSIC algorithm is designed for direction of arrival estimation in UCA and is based on the Davis Transform.

- [Davis Transforms (Davis_Transform.m) and (Davis_Transform_p1.m)]: These files represent different versions and forms of the Davis Transform implementation.

- [Dual Davis Transform (Dual_Davis_Transform.m)]: This includes various dual optimization forms of the Davis Transform implementation.

- [Auxiliary Variable Manifold Separation Technique (AV-MST) (AV_MST_v1.m), (AV_MST_V2.m), (AV_MST_V3.m)]: These implementations focus on estimating the direction of arrival using the AV-MST method. This technique models the steering vector of arbitrary array structures as the product of a sampling matrix (dependent only on the array structure) and two Vandermonde-structured wavefield coefficient vectors (dependent on the wavefield).

- [ExTG (ExTG.m)]: This version is specifically designed for estimating the optimal covariance matrix used in direction of arrival estimation.

Please explore these resources for in-depth insights into direction of arrival estimation techniques using circular array antennas. If you have any questions or require further clarification or discussion, please feel free to reach out.
