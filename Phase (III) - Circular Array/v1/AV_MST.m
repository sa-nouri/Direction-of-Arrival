%% Auxiliary Variable Manifold Seperation Technique

close all; clear all; clc

%% Matrix Covariance Calculation &  Parameter Initializing

Num_Sample = 1; SNR_dB = 100; SNR = 10 .^ (SNR_dB/10) ;
M = 11 ; d = 0: M-1 ;  Theta = [0];  Phi = [24]';
Num_Source = length(Theta); gamma = 2*pi/M * d' * 180 /pi;

ka = pi/2;
[R ,~] = CirCov(d,Theta,Phi,Num_Source,Num_Sample,SNR_dB, ka);

%% AV_MST code

K = 0:M-1; Theta = randi([0 90],length(K),1); Phi = randi([0 180],length(K),1);

% Wavefield Coefficient Vector Calculation

V_Theta = exp( 1i .* ( -M : M ) .* ( Theta * pi/180 ) ).';
V_Phi = exp( 1i .* ( -M : M ) .* ( Phi * pi/180 ) ).';
V = kron(V_Theta,V_Phi);

% Sampling Matrix Calculation
 
G_Total = [];

 for i = 1 : length(gamma)
    for j = 1 : ( 2*M + 1) 
        G(i,j) = exp(-j * (M+1-j) * gamma(i)) .* besselj(-( M+1-j ),ka);
    end
    temp = kron(G(i,:),G(i,:));
    G_Total = [G_Total; temp];

 end


% Propagator Operator and Q-Matrix Calculation

%% Music Test

angle = (0:0.1:180).';
a = exp(1i .* ka .* sind((angle' - gamma)) .* cosd(Theta));
F_Music = 1./sum(conj(a) .* (V * V') * a) ;
%% Circular Covariance Matrix
function [ R,A ] = CirCov(d,Theta,Phi,Num_Source,Num_Sample,SNR_dB, ka)
    
    M = length(d);
    gamma = 2*pi/M * d' * 180 /pi;
    A = exp(1i .* ka .* sind((Phi' - gamma)).*cosd(Theta));
    SNR = 10 .^ (SNR_dB/10) ; 
    
    noise  = (randn(length(d),Num_Sample) + ( 1i * randn(length(d),Num_Sample))) * (sqrt(1/2));
    signal = (randn(Num_Source,Num_Sample) + ( 1i * randn(Num_Source,Num_Sample))) .* (sqrt(SNR(:,1)/2));
    
    x = A * signal + noise;
    R = 1/M*(x * x');
    R = R/length(signal);
    
end