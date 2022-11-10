
%% Root_Music

SNR_dB = 20 ;Tetha = [ 66 ;88 ];
Num_Sample = 1e3; Num_Source = length(Tetha);  SNR = 100 .^ (SNR_dB/10) ; M = 10 ;d = 0: M-1; 

[R,A] = Covariance (d,Tetha,Num_Source,Num_Sample,SNR_dB);

[Betha , lambda ] = eig(R);
[~ , lambda_index ] = sort(diag(lambda),'descend');
Betha = Betha(:,lambda_index); 
U = Betha(:,Num_Source+1:end);

C = A' * U * U' * A;
Phi = sum(C,2);
Tetha_e = asin(1/pi * angle(Phi));




%%





