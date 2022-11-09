%% MVDR

SNR_dB = 20 ;Tetha = 44;
Num_Sample = 1e3; Num_Source = 1;  SNR = 100 .^ (SNR_dB/10) ; M = 10 ;d = 0: M-1; 

[R,A] = Covariance(d, Tetha, Num_Source, Num_Sample, SNR_dB);

angle = 0:0.01:180;
a = ( 1/sqrt(length(d)) ) * exp( -1i * d' * pi .* cosd(angle));
F_MVDR = (sum(conj(a) .* (inv(R) * a)))./(sum(conj(a) .* ((inv(R).^2) * a)));

plot(angle, abs(F_MVDR));

%%