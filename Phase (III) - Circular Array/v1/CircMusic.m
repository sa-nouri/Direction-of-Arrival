clc; clear; close all;

%% Parameters
Num_Sample = 1000; SNR_dB = 30; SNR = 10 .^ (SNR_dB/10) ;
M = 11 ; d = 0: M-1 ;  Theta = [10];  Phi = [10];
Num_Source = length(Theta);

%% Covariance Calculation
ka = pi;
[R ,A] = CirCov(d,Theta,Phi,Num_Source,Num_Sample,SNR_dB, ka);
for i = 1: ceil((M-1)/2)
    x(i) = abs( besselj(i-M,ka)/besselj(i,ka) );
end


%% MUSIC
[beta , lambda] = eig(R);
[~ , lambda_index ] = sort(diag(lambda),'descend');
beta = beta(:,lambda_index);
U = beta(:,Num_Source+1:end);
gamma_angle = 2*pi/M * d' * 180 /pi;
theta_vec = 0:1:90;
phi_vec = -180:1:180;
F_Music = zeros(numel(phi_vec),numel(theta_vec));
for i = 1:length(phi_vec)
    for j=1:length(theta_vec)
        a = ( 1/sqrt(length(d)) ) * exp( -1i .* ka .* sind(phi_vec(i) - gamma_angle)*cosd(theta_vec(j)));
        F_Music(i,j) =(1./sum(conj(a) .* (U * U' * a)));
    end
end
surf(theta_vec,phi_vec, 10* log10(abs(F_Music)));

%% Circular Covariance Matrix
function [ R,A ] = CirCov(d,Tetha,Phi,Num_Source,Num_Sample,SNR_dB, ka)
    
    M = length(d);
    gamma = 2*pi/M * d' * 180 /pi;
    A = exp(1i .* ka .* sind((Phi' - gamma)).*cosd(Tetha)  );
    SNR = 10 .^ (SNR_dB/10) ; 
    
    noise  = (randn(length(d),Num_Sample) + ( 1i * randn(length(d),Num_Sample))) * (sqrt(1/2));
    signal = (randn(Num_Source,Num_Sample) + ( 1i * randn(Num_Source,Num_Sample))) .* (sqrt(SNR(:,1)/2));
    
    x = A * signal + noise;
    R = 1/M*(x * x');
    R = R/length(signal);
    
end