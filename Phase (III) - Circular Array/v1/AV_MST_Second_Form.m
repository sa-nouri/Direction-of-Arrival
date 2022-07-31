%% Auxiliary Variable Manifold Seperation Technique

close all; clear; clc;

%% Matrix Covariance Calculation &  Parameter Initializing

Num_Sample = 1; SNR_dB = 100;  SNR = 10 .^ (SNR_dB/10) ; 
L = 11;     % Number of Antenna
d = 0:L-1;  Theta = [0];  Phi = [20]'; Num_Source = length(Theta); 
gamma = 2*pi/L * d' * 180 /pi;

ka = pi/2;
[R , A] = CirCov(d,Theta,Phi,Num_Source,Num_Sample,SNR_dB, ka);

%% AV_MST code

K = Num_Source; M = 100;

% Wavefield Coefficient Vector Calculation

% V_Theta = exp( 1i .* ((M-1)/2:-1:(-M+1)/2) .* ( Theta * pi/180 ) ).';
% V_Phi = exp( 1i .* ((M-1)/2:-1:(-M+1)/2) .* ( Phi * pi/180 ) ).';
% V = kron(V_Theta,V_Phi);

% Sampling Matrix Calculation

G_Total = [];
j = 0;
 for i = 1:length(gamma)
    for j = 1:M
        G(i,j) = exp(-1i * ((-M+1)/2-1+j) * gamma(i)*pi/180) .* besselj(((-M+1)/2-1+j),ka);
    end
    temp = kron(G(i,:), (G(i,:)));
    G_Total = [G_Total; temp];

 end


% Propagator Operator and Q-Matrix Calculation


%% Second Form of Sampling Matrix
% WAVEFIELD MANIFOLD SEPARATION TECHNIQUE

M = 121; Qe = 179; Qa = 90;

Theta_cal = linspace(0,90,Qa);
Phi_cal = linspace(-180,180,Qe).';

Antenna_cell = {}; V_cell = {}; Antennal_cell = {};

for k = 1 : length(gamma)
%             V_Theta = exp( 1i .* ((M-1)/2:-1:(-M+1)/2) .* ( Theta_cal(j) * pi/180 ) ).';
%             V_Phi = exp( 1i .* ((M-1)/2:-1:(-M+1)/2) .* ( Phi_cal(i) * pi/180 ) ).';
%             temp = kron(V_Phi,V_Theta);
%             V = [V, temp];
    A_temp = exp(1i .* ka .* sind((Phi_cal - gamma(k))).*cosd(Theta_cal));
    Al_temp = exp(1i .* ka .* sind(((-1*Phi_cal) - gamma(k))).*cosd(Theta_cal + 180));
    Al_temp(1,:) = []; Al_temp(end,:) = [];
    Antenna_cell{end+1} = A_temp;
    Antennal_cell{end+1} = Al_temp;
%             V_cell{end+1} = V;
end
%%

ANT_CELL = {};
G_Tot = [];

for i = 1:length(Antenna_cell)
    ANT_CELL{end+1} = [Antenna_cell{1,i}; Antennal_cell{1,i}];
   
    G = fftshift(ifft2(ANT_CELL{1,i}));
    G = G(:).';
    G_Tot = [G_Tot; G];
end

G_Total = G_Tot(:,((length(G_Tot)+1)/2 - (( M^2 -1 )/2) : (length(G_Tot)+1)/2  + (( M^2 -1)/2) ) );

%% Music Test

[beta , lambda] = eig(R);
[~ , lambda_index ] = sort(diag(lambda),'descend');
beta = beta(:,lambda_index);
U = beta(:,Num_Source+1:end);

theta_vec = linspace(0,90, 90);
phi_vec = linspace(-180, 180, 180);
% theta_vec = theta_vec0 + phi_vec0;
% phi_vec = theta_vec0 - phi_vec0;

F_Music = zeros(numel(phi_vec),numel(theta_vec));

for i = 1:length(phi_vec)
    for j=1:length(theta_vec)
        V_Theta = exp( 1i .* ((M-1)/2:-1:(-M+1)/2) .* ( theta_vec(j) * pi/180 ) ).';
        V_Phi = exp( 1i .* ((M-1)/2:-1:(-M+1)/2) .* ( phi_vec(i) * pi/180 ) ).';
        V = kron(V_Phi,V_Theta);
        
        a = G_Total*V;
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