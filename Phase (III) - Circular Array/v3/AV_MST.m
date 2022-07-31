%% Auxiliary Variable Manifold Seperation Technique

close all; clear; clc; 

%% Matrix Covariance Calculation &  Parameter Initializing

Num_Sample = 1; SNR_dB = 1000; SNR = 10 .^ (SNR_dB/10) ;
L = 11 ; d = 0:L-1 ;  Theta = [0];  Phi = [45]';

Num_Source = length(Theta); gamma = 2*pi/L * d' * 180 /pi;

ka = pi/2;
[R , A] = CirCov(d,Theta,Phi,Num_Source,Num_Sample,SNR_dB, ka);

%% AV_MST code

% K = 0:M-1; Theta = randi([0 90],length(K),1); Phi = randi([0 180],length(K),1);
K = Num_Source;

% Wavefield Coefficient Vector Calculation
M = 101;
% V_Theta = exp( 1i .* ((M-1)/2:-1:(-M+1)/2) .* ( Theta * pi/180 ) ).';
% V_Phi = exp( 1i .* ((M-1)/2:-1:(-M+1)/2) .* ( Phi * pi/180 ) ).';
% V = kron(V_Theta,V_Phi);

% Sampling Matrix Calculation
 
G_Total = [];
j = 0;
 for i = 1 : length(gamma)
    for j = 1:M
        G(i,j) = exp(-1i * ((-M+1)/2-1+j) * gamma(i)*pi/180) .* besselj(((-M+1)/2-1+j),ka);
    end
    temp = kron(G(i,:), (G(i,:)));
    G_Total = [G_Total; temp];

 end


% Propagator Operator and Q-Matrix Calculation

%% Music Test

[beta , lambda] = eig(R);
[~ , lambda_index ] = sort(diag(lambda),'descend');
beta = beta(:,lambda_index);
U = beta(:,Num_Source+1:end);
% gamma_angle = 2*pi/L * d' * 180 /pi;
theta_vec = 0:1:90;
phi_vec = -180:1:180;
F_Music = zeros(numel(phi_vec),numel(theta_vec));

for i = 1:length(phi_vec)
    for j=1:length(theta_vec)
        V_Theta = exp( 1i .* ((-M+1)/2:1:(M-1)/2) .* ( theta_vec(j) * pi/180 ) );
        V_Phi = exp( 1i .* ((-M+1)/2:1:(M-1)/2) .* ( phi_vec(i) * pi/180 ) );
        V = kron(V_Phi,V_Theta);
        a = G_Total*V.';
        F_Music(i,j) =(1./sum(conj(a) .* (U * U' * a)));
    end
end

surf(theta_vec,phi_vec, 10* log10(abs(F_Music)));


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