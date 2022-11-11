## Circular Array Antenna

The Unifrom Circular Array (UCA) has it’s elements uniformly placed along a circle of radius R. Placing the origin of the coordinate system at the center of the array. Since also this array may suffer from grating lobes, it is customary to choose the radius such that the minimum element separation is close to λ/2, as in the case of a ULA.
The following lines give particular information about each implemented matlab file for circular array antenna.

- [Ciruclar-MUSIC](CircMusicDavis.m): It is the two-dimensional multiple signal classification (MUSIC) method for UCA.

- [Circular-Vandermond](Circualr_vandermond.m): The vandermond decomposition apporach for estimating the direction of arrival for UCA.

- [Circular-MUSIC-Davis](CircMusicDavis.m): The version of MUSIC implementation for direction of arrival estimation for UCA based on Davis Transform.

- [Davis-Transforms](Davis_Transform.m),[p1](Davis_Transform_p1.m) are different versions and forms of the Davis Transform implementation.

- [Dual-Davis-Transform](Dual_Davis_Transform.m): It is various dual optimization forms of Davis Transform Implementation.

-[AUX-Var-Man-Sep_Tech](AV_MST_v1.m)[v2](AV_MST_V2.m)[v3](AV_MST_v3.m): Its primary focus is to implement the Auxiliary Variable Manifold Seperation Technique (AV-MST) for estimating direction of arrival. It is used to model the steering vector of arbitrary array structure as the product of a sampling matrix (dependent only on the array structure) and two Vandermonde-structured wavefield coefficient vectors (dependent on the wavefield).

- [ExTG](ExTG.m): It is a particular version of estimating the optimum covariance matrix to estimate the direction of arrival.
