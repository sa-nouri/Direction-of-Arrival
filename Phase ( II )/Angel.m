%%
clc ; clear ; close all

%% Angle Estimation 

Num_source = 4 ;
Num_sample = 1e3; 
SNR = [1e4 1e4 1e4 1e4 ];

M = 10 ;  % Num of Antennas 
d =(0:M-1) ; 
theta = [23 ; 40 ; 60 ; 90 ] ;%rand(Num_source,1) .* 180 / pi ;

noise  = (randn(M,Num_sample) + ( 1i * randn(M,Num_sample))) * (sqrt(1/2));
signal = (randn(Num_source,Num_sample) + ( 1i * randn(Num_source,Num_sample))) .* (sqrt(SNR(:,1)/2));

x = ( 1/sqrt(M) ) * exp( -1i .* d' .* pi .* cosd(theta)') * signal + noise;
R = x * x';
R = R/length(signal);

angle = 0:0.01:180;

a = ( 1/sqrt(M) ) * exp( -1i * d' * pi .* cosd(angle));
f_beam = diag( a' * R * a);
f_capon = diag(1./(a' * inv(R) * a));

%% finding peaks

[Pks_beam,loc_beam ] = findpeaks(abs(f_beam));
[Pks_capon,loc_capon ] = findpeaks(abs(f_capon));

angle_beam = angle(loc_beam) ;
angle_capon = angle(loc_capon) ;

%% Angel Error Calculation

error_beam = theta - angle_beam ;
error_capon = theta - angle_capon ; 

%% RMSE

rmse_beam = sqrt(diag(error_beam).^2); 
rmse_capon = sqrt(diag(error_capon).^2);

plot((rmse_beam),'b');
grid on ;
hold on ;
plot((rmse_capon),'r');
title(' RMSE  ' ) ;

%% rmse 







