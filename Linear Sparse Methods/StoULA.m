%% 
close all; clear ; clc;

%% Parameter Initializing

SNR_dB = 20 ;Theta = [ 40 ; 60 ] ;
Num_Sample = 1; Num_Source = length(Theta) ; M = 10; d = 0: M-1; 
[R, A] = Covariance(d,Theta,Num_Source,Num_Sample,SNR_dB);

%% Interpolated Array 

d1 = [ 0 2 4 6 8 9 ];
d = 0: M-1;

[R1, As] = Covariance(d1, Theta, Num_Source, Num_Sample, SNR_dB);

% As = exp( -1i .* d1' .* pi .* cosd(Theta)');
Av = interp1(d1, As, d) ; % use interp for d ---> to ULA

V = Av * As' * (inv(As * As')); % Using Text-book formoola
Rv = V * R1 * V'; % New Covariance Matrix with ULA

%% Spectrum Function

% Calculatin Music Parameters to compute Music Spectral at Text-book
[Betha , lambda ] = eig(Rv);
[~ , lambda_index ] = sort(diag(lambda), 'descend');
Betha = Betha(:, lambda_index); 
U = Betha(:, Num_Source + 1: end);
angle = 0:1:180 ;
a = (1/sqrt(length(d))) * exp( -1i * d' * pi .* cosd(angle));
    
%% Music form

% F_v = a' * 1./angle  * V' * U' * U * V * a ; % Error
Fv_music = (1./sum(conj(a) .* (U * U' * a)));

%% CVX Model

% Can not use to Multiple Sanpshot
cvx_begin
    variable Vv(length(d), length(d1)) complex
    cvx_solver sedumi
    minimize (norm((Av - (V * As)), 'fro'));
cvx_end


% *** I think output of 2 part is not same *** ?
%%  Manifold Seperation


% I think we use this formula for M (Num of Antenna) in this scope
% Other steps are same with Interpolated
% But I don't realize this scope currectly ...

M_ms = ceil( 8 * pi * max(d) / 2) ;
dd = 0 : M_ms -1 ; 

A_ms = interp1(d1,As,dd) ;
V = A_ms * As' * (inv(As * As')) ;
Rv_ms = V * R1 * V' ;    

[Betha , lambda ] = eig(Rv_ms);
[~ , lambda_index ] = sort(diag(lambda),'descend');
Betha = Betha(:,lambda_index); 
U = Betha(:, Num_Source + 1: end);
angle = 0:1:180;
aa = ( 1/sqrt(length(d)) ) * exp( -1i * dd' * pi .* cosd(angle));

Fv_ms = (1./sum(conj(aa) .* (U * U' * aa)));
    
%% Gridless
% initializing --- Data Model

Theta = [ 40 ; 60 ];    SNR_dB = 20;    Num_Sample = 1; M = 10;
SNR = 10 .^ (SNR_dB/10);    d = 0: M-1; Num_Source = length(Theta);

noise  = (randn(length(d), Num_Sample) + (1i * randn(length(d), Num_Sample))) * (sqrt(1/2));
signal = (randn(Num_Source, Num_Sample) + (1i * randn(Num_Source, Num_Sample))) .* (sqrt(SNR(:, 1)/2));
z = (1/sqrt(length(d))) * exp(-1i .* d' .* pi .* cosd(Theta)')* signal;

y = (1/sqrt(length(d))) * exp(-1i .* d' .* pi .* cosd(Theta)') * signal + noise;
Rr = y * y';
Rr = Rr/length(signal);

% Problem :: If use one snapshot, covariance matrix is not PSD.
% so we have gotten problem in next part ---- Command : chol(Rr)
% finally , we use multiple snapshot in this part to have PSD Covariance
% Matrix ---- then use one snapshot for other calculating 


%% Data Model

V = chol(Rr);
[~ , lambdaa ] = eig(V);

V = chol(Rr - (min(lambdaa) * eye(length(Rr))));
% [Betha , lambdaa ] = eig(V);

Vn = V;    V1 = V;
Vn(:, end) = []; V1(:,1) = [];

[Q, value_eig] = eig(Vn' * V1 , Vn' * Vn);

Theta = acos(1/pi * imag(log(diag(value_eig))));


%%  Determenistic Method

lambda = 25;    Etha = 1;% Reqularization Parameter and Noise Variance 

% Form (1) 

cvx_begin 
    variable z(length(d), Num_Sample) complex
    minimize (sum(norms(z, 2, 2)));
    subject to 
        norms((z - y), 2) <= Etha;
cvx_end

% Form (2) 

cvx_begin
    variable zz(length(d), Num_Sample) complex
    minimize ((lambda * sum(norms(zz, 2, 2))) + 0.5 * norms((z - y), 2));
cvx_end


%% Atomic norm
lambda = 25 ;
        
cvx_begin sdp
    variable u(1,length(d)) complex
    variable x 
    variable Z(length(d),Num_Sample) complex
    variable T(length(u),length(u)) complex
    T == toeplitz(u) ;
    cvx_solver sdpt3
    minimize(lambda/2 *( x + real(trace(T))) + 0.5 * norms((Z - y), 2));
    subject to 
        [x  Z';...
            Z  T] >= 0 ;
cvx_end

%% Henkel based

Theta = [40 ; 60];  SNR_dB = 20;    Num_Sample = 1; M = 10;
SNR = 10 .^ (SNR_dB/10);   d = 0: M-1;   Num_Source = length(Theta);

[R, A] = Covariance(d, Theta, Num_Source, Num_Sample, SNR_dB);

% m + n = length(d)  + 1 ;

m = 6;  n = 5;   dm = 0: m-1;  dn = 0: n-1;

cvx_begin 
    variable H(m, n) complex 
    M = exp(1i .* pi .* dm' .* cosd(Theta'));
    N = exp(1i .* pi .* dn' .* cosd(Theta'))';
    H =  hankel(M,N);
    variable Q1(m, m) complex
    variable Q2(Num_Sample, Num_Sample) complex
    cvx_solver sdpt3
    minimize(1/2 * real((trace(Q1) + trace(Q2))));
    subject to 
        [Q1 H';...
          H Q2] >= 0;
cvx_end


%%% *********************************************************************%%

%%
% In Hankel Based Part --- I think that I have some problem ::
% 1. m , n are cvx variables ?
% 2. How to Multiple Signal to hankel Matrix {H(M,N)} 
% 3. Is coding true ?
% 4. I use these papers for solving this part: 
% 4.1: Guaranteed minimum-rank solutions of linear matrix equations via nuclear norm minimization
% for Nuclear Norm 
% 4.2: Y. Chen, Y. Chi, Robust spectral compressed sensing via structured matrix completion,
% IEEE Trans. Inf. Theory 60 (10) (2014) 6576�6601.
