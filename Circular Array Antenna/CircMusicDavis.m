clc; clear; close all;

%% Parameters
Num_Sample = 2000; SNR_dB = 25; SNR = 10 .^ (SNR_dB/10) ;
M = 11 ; d = 0: M-1 ;  Theta = [0];  Phi = [30];
Num_Source = length(Theta);

%% Covariance Calculation
ka = pi/2;
[R ,A] = CirCov(d, Theta, Phi, Num_Source, Num_Sample, SNR_dB, ka);

for i = 1: ceil((M-1)/2)
    x(i) = abs( besselj(i-M, ka)/besselj(i, ka) );
end
alpha = 0.02;
[~,h] = max(x .* (x<=alpha));
w = exp(1i * 2 * pi /M);
F = zeros(2*h+1, M);
for i = 1:  (2*h) + 1
    for j = 1: M
        F(i,j) = w ^ ((i-h-1) * (j-1));
    end
end

J = zeros(2*h+1);
for i = 1:length(J)
    J(i, i) = 1/(1i^(i-h-1)*besselj(i-h-1, ka) * sqrt(M));
end

T = J * F ;
Rr = T * R * T';

%% MUSIC
[beta , lambda] = eig(Rr);
[~ , lambda_index ] = sort(diag(lambda),'descend');
beta = beta(:, lambda_index);
U = beta(:, Num_Source + 1: end);
gamma_angle = 2*pi/M * d' * 180 /pi;
phi_vec = -90:1:90;

% a = ( 1/sqrt(length(d)) ) * exp(1i * [-h:h]'* (phi_vec*pi/180));
% F_Music =(1./sum(conj(a) .* (U * U' * a)));

F_Music = zeros(size(phi_vec));
for i = 1:numel(phi_vec)
    a = (1/sqrt(2*h+1) ) * exp(1i * [-h:h]'* (phi_vec(i)*pi/180));
    F_Music(i) = 1/abs(a' * U * U' * a);
end
plot(phi_vec, 10 * log10(F_Music));

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