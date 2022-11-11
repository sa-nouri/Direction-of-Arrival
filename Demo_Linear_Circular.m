%% Final Project

close all; clear; clc;
%% Parameter

SNR_dB = 40 ;Theta = [ 60 80 ].';
Num_Sample = 1; Num_Source = length(Theta) ; M = 10 ;d = 0: M-1; 

A = exp(-1i .* d' .* pi .* cosd(Theta'));
SNR = 10 .^ (SNR_dB/10) ;

noise = (randn(length(d), Num_Sample) + (1i * randn(length(d), Num_Sample))) * (sqrt(1/2));
signal = (randn(Num_Source, Num_Sample) + (1i * randn(Num_Source, Num_Sample))) .* (sqrt(SNR(:,1)/2));

x = (1/sqrt(length(d)) ) * A * signal + noise;
R = x * x';
R = R/length(signal);

%% MVDR

angle = 0:0.01:180;
a = (1/sqrt(length(d))) * exp( -1i * d' * pi .* cosd(angle));
F_MVDR = (sum(conj(a) .* (inv(R) * a)))./(sum(conj(a) .* ((inv(R).^2) * a)));

%% MVDR Plot

plot(angle, -10*log10(abs(F_MVDR)));
grid on;
xlabel('Angles in degree (DOA)');
ylabel('Power (dB)');
title('MVDR (Capon) method for DOA Estimation');

%% MUSIC

[Betha , lambda ] = eig(R);
[~, lambda_index] = sort(diag(lambda), 'descend');
Betha = Betha(:, lambda_index); 
U = Betha(:, Num_Source+1:end);
angle = 0:0.01:180 ;
a = ( 1/sqrt(length(d)) ) * exp( -1i * d' * pi .* cosd(angle));
F_Music =(1./sum(conj(a) .* (U * U' * a)));
  
%% MUSIC plot

plot(angle, 10*log10(abs(F_Music)));
grid on;
xlabel('Angles in degree (DOA)');
ylabel('Power (dB)');
title('MUSIC method for DOA Estimation using Atomic Norm');
    
%% Vandermonde Decompostion

V = chol(R) ;
lambdaa = eig(V);

Vn = V ; V1 = V ;
Vn(:,end) = [] ; V1(:,1) = [] ;

[Q , value_eig ] = eig(Vn' * V1 , Vn' * Vn) ;
[~,Idx] = sort(abs(diag(value_eig)), 'descend');
val_eig_unsort = diag(value_eig);
val_eig_sort = val_eig_unsort(Idx);

Theta = acosd(1/pi * imag(log(val_eig_sort) ) )

%% Atomic Norm
y = x;
lambda = 25;
        
cvx_begin sdp
    variable u(1, length(d)) complex
    variable x 
    variable Z(length(d), Num_Sample) complex
    variable T(length(u), length(u)) complex
    T == toeplitz(u);
    cvx_solver sdpt3
    minimize( lambda/2 *( x + real(trace(T))) + 0.5 * norms((Z - y ) , 2));
    subject to 
        [x  Z';...
            Z  T] >= 0 ;
cvx_end

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
[R , A] = CirCov(d, Theta, Phi, Num_Source, Num_Sample, SNR_dB, ka);

%% Second Form of Sampling Matrix

% WAVEFIELD MANIFOLD SEPARATION TECHNIQUE
M = 200; Qe = 200; Qa = 200;
% Theta_cal = (0:Qe-1)'/Qe*180;
% Phi_cal = (0:Qa-1)/Qa*360-180;
Theta_cal = linspace(0, 180, Qe)';
Phi_cal = linspace(0, 360, Qa);
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
G_Total = G_Tot(:, L1/2-L2/2:L1/2+L2/2-1);

%% Music Test

[beta, lambda] = eig(R);
[~, lambda_index] = sort(diag(lambda), 'descend');

beta = beta(:, lambda_index);
U = beta(:, Num_Source+1:end);

theta_vec = linspace(0, 90, Qe);
phi_vec = linspace(-180, 180, Qa);

% theta_vec = theta_vec0 + phi_vec0;
% phi_vec = theta_vec0 - phi_vec0;

F_Music = zeros(numel(phi_vec), numel(theta_vec));

for i = 1:length(phi_vec)
    for j=1:length(theta_vec)
        V_Theta = exp( 1i * ((M-1)/2:-1:(-M+1)/2) .* ( theta_vec(j) * pi/180 ) ).';
        V_Phi = exp(1i * ((M-1)/2:-1:(-M+1)/2) .* ( phi_vec(i) * pi/180 ) ).';
        V = kron(V_Theta, V_Phi);
        a = G_Total*V;
        F_Music(i,j) =(1./sum(conj(a) .* (U * U' * a)));
    end
end

%% plot 2-dimensional music
figure 
surf(theta_vec, phi_vec, 10*log10(abs(F_Music)));
title(' Seperation Manifold Technique --- Multiple Snapshots')
xlabel('Theta (dg)');
ylabel('Phi (dg)');
zlabel('Power (dB)');
grid on;


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

%% Covariance

function [R,A] = Covariance(d, Tetha, Num_Source, Num_Sample, SNR_dB)

    A = exp( -1i .* d' .* pi .* cosd(Tetha'));
    SNR = 10 .^ (SNR_dB/10) ;

    noise  = (randn(length(d), Num_Sample) + ( 1i * randn(length(d), Num_Sample))) * (sqrt(1/2));
    signal = (randn(Num_Source, Num_Sample) + ( 1i * randn(Num_Source, Num_Sample))) .* (sqrt(SNR(:,1)/2));

    x = ( 1/sqrt(length(d)) ) * A * signal + noise;
    R = x * x';
    R = R/length(signal);
end