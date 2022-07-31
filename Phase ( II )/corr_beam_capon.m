
%%   Estimatin Signal
clc; clear; close all

N_sample = 1e3;        % num of signal's and noise's sample
M = 10 ;               % Num of Antennas   
N_target = 2;         % Num of targets
theta = [50 80 ]; %rand(1,N_target) * 180/pi ; 

%%SNR = 1e4 .* rand(N_target,1);
SNR = [1, 1];
d =(0:M-1) ; 

R = 0;
for i = 1:N_sample
noise  = (randn(M,1) + ( 1i * randn(M,1))) * (sqrt(1/2));
signal = (randn(N_target,1) + ( 1i * randn(N_target,1))) .* (sqrt(SNR(:,1)/2));
      
x = ( 1/sqrt(M) ) * exp( -1i .* d' .* pi .* cosd(theta)) * signal + noise;
R = R + x * x';
% plot(abs(noise))
% pause
end
R =  R/N_sample;

%%
noise  = (randn(M,N_sample) + ( 1i * randn(M,1))) * (sqrt(1/2));
signal = (randn(N_target,N_sample) + ( 1i * randn(N_target,N_sample))) .* (sqrt(SNR(:,1)/2));
x = ( 1/sqrt(M) ) * exp( -1i .* d' .* pi .* cosd(theta)) * signal + noise;
R = x*x'/N_sample;
%%

%Problems 
% how to have main signal with it's samples
% if no , to be notice size of x and noise 

%% Beamformer

angle = 0:0.01:180;
% f_beam = zeros(size(angle));
% for i = 1:length(angle)
%     a = ( 1/sqrt(M) ) * exp( -1i * d' * pi .* cosd(angle(i)));
%     f_beam(i) = a' * R * a ;
% end


 a = ( 1/sqrt(M) ) * exp( -1i * d' * pi .* cosd(angle));
 f_beam = diag( a' * R * a);

figure;
plot(angle, abs(f_beam));
grid on ; 

%% Capon 

angle = 0:0.1:180;
% f_capon = zeros(size(angle));
% for i = 1:length(angle)
%     a = ( 1/sqrt(M) ) * exp( -1i * d' * pi .* cosd(angle(i)));
%     f_capon(i) = 1 ./ (a' * inv(R) * a );
% end
a = ( 1/sqrt(M) ) * exp( -1i * d' * pi .* cosd(angle));
f_capon = diag(1./(a' * inv(R) * a));
% figure;
hold all
plot(angle, abs(f_capon));
grid on ; 


%%


