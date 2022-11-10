%% Sparse Matrix Reconstruction ---- Covariance

clear ; close all ; clc;

%% Initializing

SNR_dB = 20 ;Theta = [-35 ; 10 ; 15 ; 40];
Num_Sample = 1e3; Num_Source = length(Theta) ; M = 10; d = 0: M-1; 

[R,A] = Covariance (d,Theta,Num_Source,Num_Sample,SNR_dB);

Theta_g = -90:0.1:90 ;
A_g = exp( -1i .* d' .* pi .* sind(Theta_g));

b = sracv(R,Num_Sample,d);
fm = Music(R,d,Num_Source);

%% Fig(1) 
figure
plot(Theta_g, log10(fm), 'r') ;
hold on;
grid on;
plot(Theta_g, log10(sum(b.^2, 2)), 'b');
xlabel('DOA');
ylabel('Power(dB)');


%% Fig --- RMSE

SNR = [ 10 ; 20 ; 30 ; 50 ; 75 ; 100 ; 125 ; 150 ; 175 ; 200 ; 250 ; 300 ; 400 ; 500];
Theta = [ 15 ; 30 ; 50 ; 89]; Num_Source = length(Theta);    M = 10;   d = 0: M-1; 
Num_Sample = 1e3;

for i = 1:length(SNR)
    
    SNR_dB = SNR(i,1) ;
    [R,A] = Covariance (d,Theta,Num_Source,Num_Sample,SNR_dB);
 
    f_sracv = sracv(R , Num_Sample ,d);
 %% finding peaks

    [Pks_sracv,loc_sracv ] = findpeaks(abs(sum(f_sracv .^2, 2)), 'NPeaks', Num_Source, 'SortStr', 'descend');
    angle = 0:0.1:180;
    angle_sracv = angle(loc_sracv) ;

    r(i,1) = rmse(sort(angle_sracv), Theta);
end

%%

