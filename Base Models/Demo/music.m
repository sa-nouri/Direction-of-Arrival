%% Start 
clear ; close all ; clc;
    
%% Music

Num_Source = 1 ; Num_Sample = 1e3; 
SNR_dB = [10]; M = 10; Tetha = [ 20 ];

[R,A] = Correlation(M, Tetha, Num_Source, Num_Sample, SNR_dB) ;

[Betha, lambda ] = eig(R);
[lambda, lambda_index] = sort(diag(lambda), 'descend');
Betha = Betha(lambda_index); 

angle = 0:0.01:180;
d = ( 0 : M - 1 );

U = zeros(M ,1);
for i = (1+Num_Source):M-Num_Source
    U(i,1) =  Betha(Num_Source+i, 1);
end
 
for i = 1 : length(angle)
    a = ( 1/sqrt(M) ) * exp( -1i * d' * pi .* cosd(angle(i)));    
    q(i) = (a' * (U * U') * a);
end
 
f = diag(1./q);
%  ff = diag(1 ./ (a' * U * U' * a)) ; 
% ff = diag ( 1 ./ (abs(a' * U ).^2)); 
 
f_capon = Capon(R, M);
f_beam  = Beamformer(R, M);
 
plot(angle, abs(f));
hold on ;
plot(angle, abs(f_beam),'r');
plot(angle, abs(f_capon),'g');

 

 

