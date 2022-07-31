
%%   Estimatin Signal

N = input('Please Enter The Number of Signal & Noise >    ');

noise = randn(1,N);
signal = randn(1,N);

%d_ = input('Please Enter distance between reference element and the second sensor');

 d = rand(1,N) ;

theta = randn(1,N)* 180/pi  ;

 v = signal .* ( exp( -1i .* d .* pi .* cos(theta) ) ) ;

 SNR_i = mean( ((abs(signal)).^2) ) / mean( ((abs(noise)).^2));
 
 SNR_o = (N^2) * (SNR_i);
 
F_theta = zeros(1,length(theta));
for i = 1 : length(theta) 
    F_theta(1,i) = exp( -1i * d(1,i) * cos(theta(1,i))) ;
end

 
 G_theta = (abs((1/N).*F_theta).^2);
 
 
wi = zeros(1,length(d)) ;
for i = 1 : length(d) 
    wi(1,i) = 1/sqrt(N) * exp(-1i* d(1,i)* cos(theta(1,i))) ;
end

wb = zeros(1,length(d));
for i = 1 : length(d)
    wb(1,i) = 1/sqrt(N) * transpose(exp(-1i * d(1,i)* cos(theta(1,i)))) ;
end

y = zeros(1,length(wb));
for i = 1 : length(y)
    y(1,i) = wb(1,i)' .* signal(1,i) ;
end

R = mean((abs(y)).^2) .* (ctranspose(wb)) .* (wb) ;






