%% ------------  # Davis Transform #  ******

close all; clear; clc ;

%% Covariance Calculation & Initializing

Num_Sample = 1; SNR_dB = 100; SNR = 10 .^ (SNR_dB/10) ;
M = 11 ; d = 0: M-1 ;  Theta = [0];  Phi = [20]';
Num_Source = length(Theta); gamma = 2*pi/M * d' * 180 /pi;

ka = pi/2;
[R ,~] = CirCov(d,Theta,Phi,Num_Source,Num_Sample,SNR_dB, ka);

%% Davis Transform ------- Basic Form

for i = 1: ceil((M-1)/2)
    x(i) = abs( besselj(i-M,ka)/besselj(i,ka) );
end

alpha = 0.02;
[~,h] = max(x .* (x<=alpha));
w = exp(1i * 2 * pi /M);
F = zeros(2*h+1,M);
J = zeros(2*h+1);

for i = 1:  (2*h) + 1
    for j = 1: M
        F(i,j) = w ^ ((i-h-1) * (j-1));
    end
end

for i = 1:length(J)
    J(i,i) = 1/(1i^(i-h-1)*besselj(i-h-1, ka)*sqrt(M)) ;
end

T = J * F ;
Rr = T * R * T';

%% Davis Transform ------- P2

M0 = (M-1)/2;
Phi_Grid = [0 : 0.1 : 30]';

epsilon = 1e-3 * ones(2*M0 - 1 ,1);

    cvx_begin
        variable TR( (2*M0)-1,(2*M0)+1) complex
        minimize(square_pos(norm(TR,'fro')))
        subject to 
            for i = 1 : length(Phi_Grid)
                A = exp(1i .* ka .* sind((Phi_Grid(i) - gamma)).*cosd(Theta));
                B = exp( 1i .* (-(M0 - 1):(M0 - 1)) .* (Phi_Grid(i) * pi/180) ).';
                (max(abs(real(TR * A - B )),abs(imag(TR * A - B )))) <= epsilon ;
            end
     cvx_end
    
 Rrr = TR * R * TR' ;
 
%% 

temp = Rrr
%% Vandermonde Test

Rrr = temp
Making_PSD_Coe = 1e-7 ;
Rrr = (Rrr) * Making_PSD_Coe ;
psd = diag(1e7 *diag(Rrr));
Rrr = (Rrr +  psd );
Rrr = (Rrr + Rrr')/2;
min_eig = min(eig(Rrr));
V = chol(Rrr) ;

lambdaa = eig(V);
% V = chol(Rr - ( min(lambdaa)*0.5 * eye(length(Rr)) ) ) ;

Vn = V ; V1 = V ;
Vn(:,end) = [] ; V1(:,1) = [] ;

[Q , value_eig ] = eig(Vn' * V1 , Vn' * Vn) ;
[~,Idx] = sort(abs(diag(value_eig)),'descend');
val_eig_unsort = diag(value_eig);
val_eig_sort = val_eig_unsort(Idx);

Theta = -( 1/pi * imag( log(val_eig_sort) ) ) * 180

%% Davis Transform ------ P8 

N = 100; delta = 1e-3; k = 1:N; Phi_chos = randi([-180 180],10,100);
init_val = rand(1,1)'; epsilon = 0; G = []; M0 = (M-1)/2;

for i = 1 : 10
    
    for j = 1: length(k)
         A = exp(1i .* ka .* sind((Phi_chos(i,j) - gamma)).*cosd(Theta));
         B = exp( 1i .* (-(M0):(M0)) .* (Phi_chos(i,j) * pi/180) ).';
         temp = ( A * init_val ) - B;
         if( temp >= epsilon )
             G(i,j)= temp;
         end
    end
    
end

Theta_ach = G(1,:);

    cvx_begin
        variable lambda(M,length(k)) complex
        variable x(N,2)
        variable Q(2*N, 2*N)
        for j = 1 : length(k)
            b = exp( 1i .* (-(M0):(M0)) .* (Theta_ach(1,j) * pi/180) ).';
            minimize( ( 1/2 * x.' * Q * x ) + ( b.' * lambda(1,:) ) )
            
            subject to 
                for i = 1 : length(k)
                    A = exp(1i .* ka .* sind((Theta_ach(1,j) - gamma)).*cosd(Theta));
                    ( (Q * x ) == -( A.' * lambda(1,:) ) )    ;
                end
        end
    cvx_end
    
    for i = 1 : 10
        for j = 1: length(k)
             A = exp(1i .* ka .* sind((G(i,j) - gamma)).*cosd(Theta));
             B = exp( 1i .* (-(M0):(M0)) .* (G(i,j) * pi/180) ).';
             temp = ( A * init_val ) - B;
             if( temp == 0)
                 e(i)= temp;
             end
        end
    end

    
    cvx_begin
        variable lambda(M,length(k)) complex
        variable x(N,2)
        variable Theta_s(length(k),length(k))
        variable Q(2*N, 2*N)
        for j = 1 : length(k)
            b = exp( 1i .* (-(M0):(M0)) .* (Theta_mon(1,j) * pi/180) ).';
            minimize( ( 1/2 * x.' * Q * x ) + ( b.' * lambda(1,:) ) )
            
            subject to 
                for i = 1 : length(k)
                    A = exp(1i .* ka .* sind((Theta_mon(1,j) - gamma)).*cosd(Theta));
                    ( (Q * x ) + ( A.' * lambda(1,:) ) ) == 0 ;
                end
        end
    cvx_end
    

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