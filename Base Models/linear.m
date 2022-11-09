%% Linear Prediction

SNR_dB = 20 ;Tetha = 44;
Num_Sample = 1e3; Num_Source = 1;  SNR = 100 .^ (SNR_dB/10) ; M = 10 ;d = 0: M-1; 

[R,A] = Covariance (d,Tetha,Num_Source,Num_Sample,SNR_dB);

aa = lpc(diag(R))';


angle = 0:0.01:180;
w = ( 1/sqrt(length(d)) ) * exp( -1i * d' * pi .* cosd(angle));

H = repmat(aa,1,length(angle)) .* w ;

F_Linear = 1./(abs(sum(H, 1).^2)) ;

plot(angle, (F_Linear));

%%