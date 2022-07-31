
%% Root_Music

SNR_dB = 20 ;Tetha = 44;
Num_Sample = 1e3; Num_Source = 1;  SNR = 100 .^ (SNR_dB/10) ; M = 10 ;d = 0: M-1; 

[R,A] = Covariance (d,Tetha,Num_Source,Num_Sample,SNR_dB);

[Betha , lambda ] = eig(R);
[~ , lambda_index ] = sort(diag(lambda),'descend');
Betha = Betha(:,lambda_index); 
U = Betha(:,Num_Source+1:end);

angle = 0:0.01:180 ;
% a = ( 1/sqrt(length(d)) ) * exp( -1i * d' * pi .* cosd(angle));

C = sum((U * U'),1);
z = C * exp( -1i * d' * pi .* cosd(angle));

t = arcsin((1./(2*pi .* d) * arg(z))) ;

%%





