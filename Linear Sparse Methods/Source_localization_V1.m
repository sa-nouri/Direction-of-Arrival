%% ***** Source localization ---- CVX Algorithm

close all; clear; clc ;

%%  Parameter Initializing

Num_Sample = 1; Num_Source = 2; SNR_dB = 50; SNR = 100 .^ (SNR_dB/10) ;
M = 10 ; d = 0: M-1 ;  Tetha = [40 ; 60];

A = exp( -1i .* d' .* pi .* cosd(Tetha'));

noise  = (randn(length(d), Num_Sample) + ( 1i * randn(length(d), Num_Sample))) * (sqrt(1/2));
signal = (randn(Num_Source, Num_Sample) + ( 1i * randn(Num_Source, Num_Sample))) .* (sqrt(SNR(:,1)/2));
y = ( 1/sqrt(length(d)) ) * A * signal + noise;

Tetha_g = -90:0.1:90 ;
A_g = exp( -1i .* d' .* pi .* sind(Tetha_g));

%% CVX --- l1-magic

%% Matlab version

x = (A_g)' * (y);

%% CVX Version

fprintf(1, 'Single Sample Source locarization with l1');

cvx_begin quiet
    
    variable x_l1(length(Tetha_g), Num_Sample) complex
    minimize(norm(x_l1, 1))
    subject to 
            A_g * x_l1 == y;
cvx_end

%% Iterative 

Max_Itr = 100 ;
W = rand(length(Tetha_g), Num_Sample); % initial Weights
lambda = 0.5 ; 

for k = 1: Max_Itr
    cvx_begin quiet
    %cxv_solver sedumi
        variable x_log(length(Tetha_g), Num_Sample) complex
        minimize ( sum(W .* abs(x_log)))
        % minimize ( 0.5 * sum_squre( y - A_g * x_log) + lambada * norm(x_log , 1) 
        subject to
            real(A_g) * real(x_log) <= real(y);
    cvx_end
end


%% CVX -- rank Reduction

% Dimensionality Reduction
ReducedDim = 4;
y_svd = mySVD(y,ReducedDim);

% CVX version

cvx_begin
    variable x_d(length(Tetha_g), Num_Sample) complex
    minimize ( lambda*(norm(x_d,2)) + 0.5 * (norm(A_g * x_d - y_svd , 'fro' )));
cvx_end

% ***************** Finish *********************** %



