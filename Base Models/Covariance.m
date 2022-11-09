%% Covariance

function [R,A] = Covariance(d, Tetha, Num_Source, Num_Sample, SNR_dB)

    A = exp( -1i .* d' .* pi .* cosd(Tetha'));
    SNR = 100 .^ (SNR_dB/10) ;

    noise  = (randn(length(d),Num_Sample) + ( 1i * randn(length(d),Num_Sample))) * (sqrt(1/2));
    signal = (randn(Num_Source,Num_Sample) + ( 1i * randn(Num_Source,Num_Sample))) .* (sqrt(SNR(:,1)/2));

    x = ( 1/sqrt(length(d)) ) * exp( -1i .* d' .* pi .* cosd(Tetha')) * signal + noise;
    R = x * x';
    R = R/length(signal);

end

%% 