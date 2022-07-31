%% Optimization Problem to find G

close all; clear; clc
%% Parameter Initializing

Num_Sample = 1; SNR_dB = 50; SNR = 10 .^ (SNR_dB/10) ;
M = 11 ; d = 0: M-1 ;  Theta = [0];  Phi = [20]';
Num_Source = length(Theta); gamma = 2*pi/M * d' * 180 /pi; ka = pi/2;

[R ,~] = CirCov(d,Theta,Phi,Num_Source,Num_Sample,SNR_dB, ka);

%% Davis Transform ------- P2 Optimization
% 
% M0 = (M-1)/2; Phi_Grid = [0 : 0.1 : 30]';
% epsilon = 1e-3 * ones(2*M0 - 1 ,1);
% 
%     cvx_begin
%         variable TR( (2*M0)-1,(2*M0)+1) complex
%         minimize(square_pos(norm(TR,'fro')))
%         subject to 
%             for i = 1 : length(Phi_Grid)
%                 A = exp(1i .* ka .* sind((Phi_Grid(i) - gamma)).*cosd(Theta));
%                 B = exp( 1i .* (-(M0 - 1):(M0 - 1)) .* (Phi_Grid(i) * pi/180) ).';
%                 (max(abs(real(TR * A - B )),abs(imag(TR * A - B )))) <= epsilon ;
%             end
%      cvx_end
%     
% Rrr = TR * R * TR' ;
% 
%     % Testing Optimization by Vandermonde Decomposition
%     [Theta, lambda] = vand_decom(Rrr);
%     
%% Manifold Seperation --- Optimization Problem + Result

    % Parameter Initializing
        T = 100; % it should be grater than length(A), and length(D)
        Phi_G = linspace(0, 360, T); M_ = (2 * M) + 1;
        
        A = exp(1i .* ka .* sind((Phi_G - gamma)) .* cosd(Theta)); % has a different with our Cov matrix
        D = exp(1i .* (-(M_ -1)/2 : (M_ + 1)/2).' .* Phi_G);        

        %         %% First Form
%         
% cvx_begin
%     variable G1(M_ + 1, M) complex
%     Objective = 0;
%     
%     for i = 1:T
%         Objective = Objective + square_pos(norm((G1*A(:,i) - D(:,i)), 'fro'));  
%     end
%     
%     minimize(Objective)
% cvx_end
% 
%     % norm(G1-G_ls) = 1.8156e-07;
%     
        %% Second Form
        
cvx_begin
        variable G2(M_ + 1, M) complex
        minimize (square_pos(norm((G2*A - D),'fro')))
cvx_end

    % norm(G2-G_ls) = 4.5885e-10;
    % Test
[Theta, lambda] = vand_decom(G2 * R * G2.')
[~, I] = max(lambda);
Theta(I)

%% Least Square Form
        
G_ls = A * D' * inv(D * D');

%% Circular Covariance Matrix
function [ R,A ] = CirCov(d,Theta,Phi,Num_Source,Num_Sample,SNR_dB, ka)
    
    M = length(d);
    gamma = 2*pi/M * d' * 180 /pi;
    A = exp(1i .* ka .* sind((Phi' - gamma)) .* cosd(Theta));
    SNR = 10 .^ (SNR_dB/10) ; 
    
    noise  = (randn(length(d),Num_Sample) + ( 1i * randn(length(d),Num_Sample))) * (sqrt(1/2));
    signal = (randn(Num_Source,Num_Sample) + ( 1i * randn(Num_Source,Num_Sample))) .* (sqrt(SNR(:,1)/2));
    
    x = A * signal + noise;
    R = 1/M*(x * x');
    R = R/length(signal);
    
end

%% Vandermonde Decomposition to Find Arrival
function [Theta, lambda] = vand_decom(Rr)

    Making_PSD_Coe = 1e-3 ; Rr = (Rr) * Making_PSD_Coe ;
    psd = diag(1e1 *diag(Rr)); Rr = (Rr +  psd );
    Rr = (Rr + Rr')/2; %min_eig = min(eig(Rr));
%     V = chol(Rr) ;

    lambda = eig(Rr);
    V = chol(Rr + ( min(lambdaa)* 0.01 * eye(length(Rr)) ) ) ;

    Vn = V ; V1 = V ;
    Vn(:,end) = [] ; V1(:,1) = [] ;

    [Q , value_eig ] = eig(Vn' * V1 , Vn' * Vn) ;
    [~,Idx] = sort(abs(diag(value_eig)),'descend');
    val_eig_unsort = diag(value_eig);
    val_eig_sort = val_eig_unsort(Idx);

    Theta = -( 1/pi * imag( log(val_eig_sort) ) ) * 180
end