%% Auxiliary Variable Manifold Seperation Technique
close all; clear; clc;
%% Matrix Covariance Calculation &  Parameter Initializing

Num_Sample = 1000;
SNR_dB = 100;
SNR = 10 .^ (SNR_dB/10) ;
L = 11;     % Number of Antenna
d = 0:L-1;
Theta = [30 60];
Phi = [10 -12]';
Num_Source = length(Theta);
gamma = 2*pi/L * d' * 180 /pi;
ka = pi/2;
[R , A] = CirCov(d,Theta,Phi,Num_Source,Num_Sample,SNR_dB, ka);

%% Second Form of Sampling Matrix

% WAVEFIELD MANIFOLD SEPARATION TECHNIQUE
M = 100; Qe = 100; Qa = 100;
% Theta_cal = (0:Qe-1)'/Qe*180;
% Phi_cal = (0:Qa-1)/Qa*360-180;
Theta_cal = linspace(0,180, Qe)';
Phi_cal = linspace(0,360, Qa);
G_Tot = [];
for idx = 1 : length(gamma)
    A_temp = exp(1i * ka * sind((Phi_cal - gamma(idx))).*cosd(Theta_cal));
    Al_temp = exp(1i * ka * sind((Phi_cal + 180 - gamma(idx))).*cosd((Theta_cal(2:end-1))));
    Al_p = [A_temp; Al_temp];
    
    G_idx = fftshift(ifft2(Al_p)).';
    G_idx = G_idx(:);
    G_Tot = [G_Tot; G_idx.'];
end

L1 = size(G_Tot,2);
L2 = M^2;
% G_Total = G_Tot(:,(length(G_Tot)+1)/2 - (( M^2 -1 )/2) : (length(G_Tot)+1)/2  + (( M^2 -1)/2));
G_Total = G_Tot(:,L1/2-L2/2:L1/2+L2/2-1);

%% Music Test

[beta , lambda] = eig(R);
[~ , lambda_index ] = sort(diag(lambda),'descend');
beta = beta(:,lambda_index);
U = beta(:,Num_Source+1:end);

theta_vec = linspace(0, 90, Qe);
phi_vec = linspace(-180, 180, Qa);
% theta_vec = theta_vec0 + phi_vec0;
% phi_vec = theta_vec0 - phi_vec0;

F_Music = zeros(numel(phi_vec),numel(theta_vec));

for i = 1:length(phi_vec)
    for j=1:length(theta_vec)
        V_Theta = exp( 1i * ((M-1)/2:-1:(-M+1)/2) .* ( theta_vec(j) * pi/180 ) ).';
        V_Phi = exp(1i * ((M-1)/2:-1:(-M+1)/2) .* ( phi_vec(i) * pi/180 ) ).';
        V = kron(V_Theta, V_Phi);
        a = G_Total*V;
        F_Music(i,j) =(1./sum(conj(a) .* (U * U' * a)));
    end
end

surf(theta_vec,phi_vec, 10*log10(abs(F_Music)));



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