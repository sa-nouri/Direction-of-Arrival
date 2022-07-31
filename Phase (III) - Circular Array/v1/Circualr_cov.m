%% 

clear ; close all  ; clc ; 

%% Parameters
Num_Sample = 1; SNR_dB = 100; SNR = 10 .^ (SNR_dB/10) ;
M = 11 ; d = 0: M-1 ;  Theta = [0];  Phi = [40];
Num_Source = length(Theta);

%% Covariance Calculation
ka = pi/2;
[R ,A] = CirCov(d,Theta,Phi,Num_Source,Num_Sample,SNR_dB, ka);

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
%% Vandermonde Test

  
epsilon = 0.1;
for i = 1:(2*h+1)
    Rr(i,i) = real(Rr(i,i));
end
min_eig = min(eig(Rr));
V = chol(Rr + epsilon*eye(2*h+1)) ;
lambdaa = eig(V);

% V = chol(Rr - ( min(lambdaa)*0.5 * eye(length(Rr)) ) ) ;
% [Betha , lambdaa ] = eig(V);

Vn = V ; V1 = V ;
Vn(:,end) = [] ; V1(:,1) = [] ;

[Q , value_eig ] = eig(Vn' * V1 , Vn' * Vn) ;
[~,Idx] = sort(abs(diag(value_eig)),'descend');
val_eig_unsort = diag(value_eig);
val_eig_sort = val_eig_unsort(Idx);

Theta = -(1/pi * imag(log(val_eig_sort) ) )*180

%% Music Test
% [Betha , lambda ] = eig(R);
% [~ , lambda_index ] = sort(diag(lambda),'descend');
% Betha = Betha(:,lambda_index); 
% U = Betha(:,Num_Source+1:end);
% gamma_angle = 2*pi/M * d' * 180 /pi;
% angle = -90:0.1:90 ;
% a = ( 1/sqrt(length(d)) ) * exp( -1i .* ka .* sind((angle - gamma_angle)));
% F_Music =(1./sum(conj(a) .* (U * U' * a)));
% plot(angle, 10* log10(abs(F_Music)));
    
% %% Davis and Atomic norm minimization
%  lambda = 25 ;
%  gamma = 0.1 ;
%         
% cvx_begin sdp
%         variable u1
%         variable u(1,length(d)-1) complex
%         variable x 
%         variable Z(length(d),1) complex
%         variable TT(M,M) complex
% %         variable T(length(u),length(u)) complex
%         T = toeplitz([u1,u]) ;
%         T = TT  * T * TT' ; 
%         cvx_solver sdpt3
%         minimize( lambda/2 *( x + (trace(T)) ) +  0.5 * norms((Z - y ) , 2)  + gamma * norms(TT,2)) ;
%         subject to 
%             [ x  Z';...
%               Z  T ] >= 0 ;
% cvx_end

%%  Cir to Rec
%% Circular Covariance Matrix

function [ R,A ] = CirCov(d,Tetha,Phi,Num_Source,Num_Sample,SNR_dB, ka)
    
    M = length(d);
    gamma = 2*pi/M * d' * 180 /pi;
    A = exp(1i .* ka .* cosd((Phi' - gamma)).*cosd(Tetha)  );
    SNR = 10 .^ (SNR_dB/10) ; 
    
    noise  = (randn(length(d),Num_Sample) + ( 1i * randn(length(d),Num_Sample))) * (sqrt(1/2));
    signal = (randn(Num_Source,Num_Sample) + ( 1i * randn(Num_Source,Num_Sample))) .* (sqrt(SNR(:,1)/2));
    
    x = A * signal + noise;
    R = 1/M* (x * x');
    R = R/length(signal);
    
end

%% 
