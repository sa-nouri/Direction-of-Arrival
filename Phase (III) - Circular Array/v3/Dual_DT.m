%% ------------  Dual Optimization for Davis Transform  ******

close all; clear; clc ;

%% Covariance Calculation & Initializing

Num_Sample = 1; SNR_dB = 100; SNR = 10 .^ (SNR_dB/10) ;
M = 11 ; d = 0: M-1 ;  Theta = [0];  Phi = [24]';
Num_Source = length(Theta); gamma = 2*pi/M * d' * 180 /pi;

ka = pi/2;
[R ,~] = CirCov(d,Theta,Phi,Num_Source,Num_Sample,SNR_dB, ka);

%% Transformation UCA in Signal Environment

% Acording to algorithm which is attached in the paper
% Coding is based on Steps 

%% Step (1) --- Initialization

N = M; x = rand(2*N,1) * 0.01;  % use inital value of paper
delta = 1e-5; I = 8; 

R = eye(N); 
Q = [real(R) -imag(R) ; imag(R) real(R)];
Epsilon = zeros(length(Q),1);
epsilon = 1e-2;

%% Itertation
i = 0;
temp = delta;

while((i < I) || (temp > delta))
    i = i +1;
    K = (1:(2^(4+i)) +1) ;
    Theta_Grid = (360)/(length(K)-1) .*( K - 1) ;
    [Theta_ach] = findG(ka,gamma, x, Theta_Grid,Epsilon,epsilon); % G is not Create.
    [X, lambda] = DQSIP8(ka,gamma, Theta_ach, Q, epsilon); % (P8k) Solution
    
    temp = abs((X.' * Q * X) - ( x.' * Q * x)) ;
    [E] = findE(ka,gamma, Theta_ach,X, epsilon);
    Epsilon = E; 
    x = X;
end

    [X, lambda,Theta_Opt] = DualQSIP7(ka,gamma,K, N, Q); % (P7k) Solution
    


%% Circular Covariance Matrix
function [ R,A ] = CirCov(d,Theta,Phi,Num_Source,Num_Sample,SNR_dB, ka)
    
    M = length(d);
    gamma = 2*pi/M * d' * 180 /pi;
    A = exp(1i .* ka .* sind((Phi' - gamma)).*cosd(Theta));
    SNR = 10 .^ (SNR_dB/10) ; 
    
    noise  = (randn(length(d),Num_Sample) + ( 1i * randn(length(d),Num_Sample))) * (sqrt(1/2));
    signal = (randn(Num_Source,Num_Sample) + ( 1i * randn(Num_Source,Num_Sample))) .* (sqrt(SNR(:,1)/2));
    
    x = A * signal + noise;
    R = 1/M*(x * x');
    R = R/length(signal);
    
end

%% Problem (6) Solution ----- Quadratic Semi-Infinite Programming
function [x] = QSIP6(ka, gamma, N, Theta_Grid, Q)
    
    epsilon = 1e-5;
    
    for i = 1: length(Theta_Grid)
        a = exp(1i .* ka .* cosd((Theta_Grid(i) - gamma))).'; % gamma value
        A = [real(a) -imag(a) ; -real(a) imag(a) ; imag(a) real(a); - imag(a) -real(a)] ;
        b = exp( 1i .* ( Theta_Grid(i) * pi/180) ).'; % "b" formation
        B = [epsilon+real(b) epsilon-real(b) epsilon+imag(b) epsilon-imag(b)].';
        cvx_begin
            variable x(2*N, 1)
            minimize(1/2 * x.' * Q * x)
            subject to
                ( A*x - B ) <= zeros(length(B),1);
        cvx_end
    end 

end

%% Problem (7) Solution ---- Dual Parametirization ---- Quadratic Semi-Infinite Programming
function [x, lambda,Theta_Opt] = DualQSIP7(ka, gamma,K, N, Q)
   
    
    cvx_begin
        variable x(2*N,1)
        variable lambda(4,length(K))
        variable Theta_Opt(length(K),1)
        b = exp( 1i .* ( Theta_Opt * pi/180) );
        B = [epsilon+real(b); epsilon-real(b); epsilon+imag(b); epsilon-imag(b)].';            
        minimize((1/2 * x.' * Q * x) + trace(B * lambda))
        subject to 
                for j = 1 : length(K)
                    a = exp(1i .* ka .* cosd((Theta_Opt(j) - gamma))).';
                    A = [real(a) -imag(a) ; -real(a) imag(a) ; imag(a) real(a); - imag(a) -real(a)] ;
                    temp = temp + ( A.' * lambda(:,j)) ;
                end
            ( Q * x ) + temp == zeros(length(Q*x),1);
        
    cvx_end

end

%% Problem (8) Solution ---- Dual Quadratic Semi-Infinite Programming
function [x, lambda] = DQSIP8(ka,gamma, Theta_Grid,Q, epsilon)
    
    k = length(Theta_Grid) ;
    temp = 0;
    cvx_begin
        variable x(2*(length(gamma)),1)
        variable lambda(4,k)

        b = exp( 1i .* ( Theta_Grid * pi/180) );
        B = [epsilon+real(b); epsilon-real(b); epsilon+imag(b); epsilon-imag(b)].';            
        minimize((1/2 * x.' * Q * x) + trace(B * lambda)); 
        
        subject to
            for j = 1 : k
                a = exp(1i .* ka .* cosd((Theta_Grid(j) - gamma))).';
                A = [real(a) -imag(a) ; -real(a) imag(a) ; imag(a) real(a); - imag(a) -real(a)].' ;
                temp = temp + sum(A * lambda) ; % sum should be trace but there is error % Size of A is not possible to have Trace
            end
            ( Q * x ) + temp == zeros(length(Q*x),1);
    
    cvx_end

end

%% Find G-Matrix 
function [G] = findG(ka,gamma, x, Theta_Grid,Epsilon,epsilon)
    
    G = Epsilon;
    for i = 1: length(Theta_Grid)
        a = exp(1i .* ka .* cosd((Theta_Grid(i) - gamma))).'; % gamma value & 
        A = [real(a) -imag(a) ; -real(a) imag(a) ; imag(a) real(a); - imag(a) -real(a)] ;
        b = exp( 1i .* ( Theta_Grid(i) * pi/180) ).'; % "b" formation
        B = [epsilon+real(b) epsilon-real(b) epsilon+imag(b) epsilon-imag(b)].';
        
        if( (( A*x - B ) >= 0)) %% There is Problem : 2 elements of matrix are Positve and others are Negative. :))
            G = [G Theta_Grid(i)];
        end 
    end
    G = unique(G);
end

%% Finda E-Matrix
function [E] = findE(ka,gamma, G,x, epsilon)
    E = [];
    for i = 1: length(G)
        a = exp(1i .* ka .* cosd((G(i) - gamma))).'; % gamma value & 
        A = [real(a) -imag(a) ; -real(a) imag(a) ; imag(a) real(a); - imag(a) -real(a)] ;
        b = exp( 1i .* ( G(i) * pi/180) ).'; % "b" formation
        B = [epsilon+real(b) epsilon-real(b) epsilon+imag(b) epsilon-imag(b)].';
        
        if( ( A*x - B ) == 0 )
            E = [ E G(i)];
        end
    end
end