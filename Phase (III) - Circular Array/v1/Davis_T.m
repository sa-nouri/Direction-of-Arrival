%% ------------
close all; clear; clc ;

%% Covariance Calculation & Initializing

Num_Sample = 1; SNR_dB = 100; SNR = 10 .^ (SNR_dB/10) ;
M = 11 ; d = 0: M-1 ;  Theta = [0];  Phi = [40];
Num_Source = length(Theta);

ka = pi/2;
[R ,A] = CirCov(d,Theta,Phi,Num_Source,Num_Sample,SNR_dB, ka);

%% Davis Transform ------- Basic Form
for i = 1: ceil((M-1)/2)
    x(i) = abs( besselj(i-M,ka)/besselj(i,ka) );
end
alpha = 0.02;
[~,h] = max(x .* (x<=alpha));
w = exp(1i * 2 * pi /M);
F = zeros(2*h+1,M);
for i = 1:  (2*h) + 1
    for j = 1: M
        F(i,j) = w ^ ((i-h-1) * (j-1));
    end
end

J = zeros(2*h+1);
for i = 1:length(J)
    J(i,i) = 1/(1i^(i-h-1)*besselj(i-h-1, ka)*sqrt(M)) ;
end
T = J * F ;
Rr = T * R * T';

%% Davis Transform ---- Fsemnif

Phi_Grid = [-180 : 0.01 : 180]';
M0 = (M-1)/2;
ntheta = (2);

objfun = @(x) (norm(x,'fro')).^2 ;
A = exp(1i .* ka .* sind((Phi_Grid' - gamma)) .* cosd(Theta));
B = exp( 1i .* (-M0:M0) .* (Phi*pi/180) ).';

x0 = 0.2 .* ( ones(2 * M0 + 1) + 1i * ones(2 * M0 + 1)) ;
x_real = fseminf(objfun,real(x0),ntheta,@seminfcon);
x_imag = fseminf(objfun,imag(x0),ntheta,@seminfcon);

x = x_real + x_imag;
    

%% Seminfcon

function [c, ceq, K1, K2, s] = seminfcon(X,s)

% Initial sampling interval
if isnan(s(1,1))
   s = rand(2,2);
end

M = length(X); d = 0 : M -1 ; ka = pi/2;
gamma = 2*pi/M * d' * 180 /pi ;
M0 = (M-1)/2;  x = [-180:0.01:180];
B = exp( 1i .* (-M0:M0) .* x' * pi/180  ).';

W1 = [-180:0.01:180]' .* s(1,1);
W2 = [-180:0.01:180]' .* s(2,1);

% Sampling set
W11 = exp(1i .* ka .* sind((W1' - gamma)) );
W22 = exp(1i .* ka .* sind((W2' - gamma)) );

Y1 = X * W11;

% Semi-infinite constraint 
K1 = norm( ((X * W11) - B),2);
K2 = norm( ((X * W22) - B),2);

% No finite nonlinear constraints
c = []; ceq=[];
end


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