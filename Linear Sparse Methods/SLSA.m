%% ***** Source localization ---- CVX Algorithm

clear ; close all ; clc;

%%  Parameter Initializing

Num_Sample = 200; SNR_dB = 10; SNR = 10 .^ (SNR_dB/10) ;
M = 8 ; d = 0: M-1 ;  Theta = [40; 80; 110; 140];
Num_Source = length(Theta);
A = exp( -1i .* d' .* pi .* cosd(Theta'));

noise  = (randn(length(d),Num_Sample) + ( 1i * randn(length(d),Num_Sample))) * (sqrt(1/2));
signal = (randn(Num_Source,Num_Sample) + ( 1i * randn(Num_Source,Num_Sample))) .* (sqrt(SNR(:,1)/2));
y = ( 1/sqrt(length(d)) ) * A * signal + noise;

Tetha_g = 0:1:180 ;
A_g = exp( -1i .* d' .* pi .* cosd(Tetha_g));


%% CVX --- l1 ---- BPDN

fprintf(1, 'Single Sample Source locarization with BPDN \n');

cvx_begin
    
    variable x_bpdn(length(Tetha_g),Num_Sample) complex
    minimize(norm(x_bpdn,1))
    subject to 
             norm((A_g * x_bpdn) - y) <= 24  ; 
cvx_end

figure()
plot(Tetha_g , 20 * log10(abs(x_bpdn)))

%% l1-norm --- Multiple snapshots

fprintf(1, 'L1-norm multiple snapshots source locarization \n');

lambda = 24 ;

    cvx_begin
        variable S(length(Tetha_g),Num_Sample)
        cvx_solver sedumi
        minimize (0.5 * square_pos(norm((A_g * S)- y , 'fro' ))+ lambda * sum(norms(S, 2, 2))); %% norms(S,P,Dim);
    cvx_end

figure()
plot(Tetha_g , 20*log10(sqrt(sum(abs(S).^2,2))))

%% CVX --- l1-SVD

fprintf(1, 'l2,1 norm optimization for multiple snapshots\n');

Y = sv(y,Num_Source,Num_Sample);
lambda = 24 ;

    cvx_begin
        variable S_sv(length(Tetha_g),Num_Source)
%         minimize ( power(2,norm((Y - A_g * S_sv ) ,'fro')) + lambda * sum(norms(S_sv)));
        minimize (0.5*square_pos(norm((A_g * S_sv)- Y , 'fro' )) + lambda * sum(norms(S_sv, 2, 2)));
    cvx_end

figure()
plot(Tetha_g , 20*log10(sqrt(sum(abs(S_sv) .^ 2, 2))))   

%%  Fig(7) -------------   SNR = 10 dB ---- M-1 Source

[R,A] = Covariance(d,Theta, Num_Source, Num_Sample, SNR_dB) ;
F_beam = Beamformer(R, d, Tetha_g);
F_capon = Capon(R, d, Tetha_g);
F_Music = Music(R, d, Num_Source) ;
F_svd = norm_lij(S_sv, 2, 1);


figure()
plot(Tetha_g,log10(abs(F_beam)),'r');
hold on ;
grid on ;
plot(Tetha_g,log10(abs(F_capon)),'b');
% plot(Tetha_g,log10(abs(F_Music)),'g');
plot(Tetha_g,log10(F_svd),'y');
xlabel('DOA (degrees)');
ylabel('Power (dB) ');
title(' 7 Sources & SNR = 10 dB' ) ;
hold off ;


%% Figure(8) --- Music & L1-svd Sensitivity to Num of Sources  

% 
% music 

F1 = Music(R,d,1) ;
F2 = Music(R,d,2) ;
F3 = Music(R,d,3) ;
F4 = Music(R,d,4) ;
F7 = Music(R,d,7) ;
% 
figure()
plot(Tetha_g,log10(abs(F1)),'g');
hold on; 
grid on;
plot(Tetha_g, log10(abs(F2)),'b');
plot(Tetha_g, log10(abs(F3)),'y');
plot(Tetha_g, log10(abs(F4)),'r');
plot(Tetha_g, log10(abs(F7)),'k');
hold off ;

xlabel('DOA (degrees)');
ylabel('Music Power (dB)');
title(' Music sensitivity for Num of Source');

% l1-svd 

x1 =  l1_svd(y,A_g , 1 ,Num_Sample);
x2 =  l1_svd(y,A_g , 2 ,Num_Sample);
x3 =  l1_svd(y,A_g , 3 ,Num_Sample);
x4 =  l1_svd(y,A_g , 4 ,Num_Sample);
x7 =  l1_svd(y,A_g , 7 ,Num_Sample);

F_svd1 = norm_lij(x1, 2, 1);
F_svd2 = norm_lij(x2, 2, 1);
F_svd3 = norm_lij(x3, 2, 1);
F_svd4 = norm_lij(x4, 2, 1);
F_svd7 = norm_lij(x7, 2, 1);

figure(9)
plot(Tetha_g,log10(F_svd1), 'y');
hold on;
grid on;
plot(Tetha_g,log10(F_svd2), 'r');
plot(Tetha_g,log10(F_svd3), 'b');
plot(Tetha_g,log10(F_svd4), 'g');
plot(Tetha_g,log10(F_svd7), 'k');

xlabel('DOA (degrees)');
ylabel('Power (dB) ');
title('l1-svd sensitivity to num of sources');
hold off;