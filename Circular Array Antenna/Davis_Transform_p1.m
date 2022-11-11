%% Davis Transform

Num_Sample = 1000; SNR_dB = 30; SNR = 10 .^ (SNR_dB/10) ;
M = 11 ; d = 0: M-1 ;  Theta = [10];  Phi = [10];
Num_Source = length(Theta);
epsilon = 1e-3;

ka = pi;
[R, A] = CirCov(d, Theta, Phi, Num_Source, Num_Sample, SNR_dB, ka);

cvx_begin
    variable TR(length(d),length(d)) complex
    B = TR * A;
    minimize(square(norm(TR,'fro')))
    subject to 
            abs(TR * A - B ) < epsilon;
cvx_end


%% Circular Covariance Matrix
function [R, A] = CirCov(d, Tetha, Phi, Num_Source, Num_Sample, SNR_dB, ka)
    
    M = length(d);
    gamma = 2 * pi/M * d' * 180 /pi;
    A = exp(1i .* ka .* cosd((Phi' - gamma)).*cosd(Tetha)  );
    SNR = 10 .^ (SNR_dB/10); 
    
    noise  = (randn(length(d), Num_Sample) + (1i * randn(length(d), Num_Sample))) * (sqrt(1/2));
    signal = (randn(Num_Source, Num_Sample) + (1i * randn(Num_Source, Num_Sample))) .* (sqrt(SNR(:, 1)/2));
    
    x = A * signal + noise;
    R = 1/M * (x * x');
    R = R/length(signal);
end