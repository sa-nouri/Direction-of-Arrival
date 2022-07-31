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

N = 11; x = rand(2*N,1);  % use inital value of paper
delta = 1e-5; I = 8; K = 1:N ;
Theta_Grid = (180)/(N-1) .*( K - 1) ;

R = eye(N, N) * ( 1 + 1i); % ?????
Q = [real(R) -imag(R) ; imag(R) real(R)];
Epsilon = zeros(length(Q),1);

%% Itertation

for i = 1:N % max of i ?
    [Theta_ach] = findG(ka,gamma, x, Theta_Grid); % G is not Create.
    [X, lambda] = DQSIP8(gamma, N, Theta_Grid,Q); % (P8k) Solution
    
    temp = abs((X(end,:).' * Q * X(end,:)) - ( x.' * Q * x)) ;
    
    if( (i < I) || (temp > delta))
        [E] = findE(ka,gamma, G,x);
    else
        [X, lambda,Theta_Opt] = DualQSIP7(ka,gamma, N, Q); % (P7k) Solution
    end
    
    x =  X(end,:);
    Thete_Grid = Theta_Opt;
end

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
function [X] = QSIP6(ka, gamma, N, Theta_Grid, Q)
    
    epsilon = 1e-5;
    X = zeros(length(Theta_Grid), 2*N);
    
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
        X(i,:) = x;
    end 

end

%% Problem (7) Solution ---- Dual Parametirization ---- Quadratic Semi-Infinite Programming
function [X, lambda,Theta_Opt] = DualQSIP7(ka, gamma, N, Q)
    
    i = 10;
    k = (2^(i+4)) +1 ;
    X = zeros(k, 2*N);
    
    cvx_begin
        variable x(2*N,1)
        variable lambda(4,k)
        variable Theta_Opt(k,1)
        value = 0;
        for j = 1 : k
            b = exp( 1i .* ( Theta_Opt(j) * pi/180) ).';
            B = [epsilon+real(b) epsilon-real(b) epsilon+imag(b) epsilon-imag(b)].';
            value = B.' * lambda(:,j);
        end
            minimize((1/2 * x.' * Q * x) + value) 
            
        temp = length(2*N,1) ;
        for j = 1 : k
            a = exp(1i .* ka .* cosd((Theta_Opt(j) - gamma))).';
            A = [real(a) -imag(a) ; -real(a) imag(a) ; imag(a) real(a); - imag(a) -real(a)] ;
            temp = temp + ( A.' * lambda(:,j)) ;
        end
        subject to 
            ( Q * x ) + temp == zeros(length(Q*x),1);
        
        X(j,:) = x;
    cvx_end

end

%% Problem (8) Solution ---- Dual Quadratic Semi-Infinite Programming
function [X, lambda] = DQSIP8(ka,gamma, N, Theta_Grid,Q)
    
    k = length(Theta_Grid) ;
    X = zeros(length(Theta_Grid), 2*N);
    
    cvx_begin
        variable x(2*N,1)
        variable lambda(4,k)
        for j = 1 : k
            b = exp( 1i .* ( Theta_Grid(j) * pi/180) );
            B = [epsilon+real(b) epsilon-real(b) epsilon+imag(b) epsilon-imag(b)];
            value = B.' * lambda(:,j);
        end
            minimize((1/2 * x.' * Q * x) + value) 
            
        temp = length(2*N,1) ;
        for j = 1 : k
            a = exp(1i .* ka .* cosd((Theta_Grid(j) - gamma))).';
            A = [real(a) -imag(a) ; -real(a) imag(a) ; imag(a) real(a); - imag(a) -real(a)] ;
            temp = temp + ( A.' * lambda(:,j)) ;
        end
        subject to 
            ( Q * x ) + temp == zeros(length(Q*x),1);
    
    cvx_end

end

%% Find G-Matrix 
function [G] = findG(ka,gamma, x, Theta_Grid)
    
    G = []; epsilon = 1e-4;
    for i = 1: length(Theta_Grid)
        a = exp(1i .* ka .* cosd((Theta_Grid(i) - gamma))).'; % gamma value & 
        A = [real(a) -imag(a) ; -real(a) imag(a) ; imag(a) real(a); - imag(a) -real(a)] ;
        b = exp( 1i .* ( Theta_Grid(i) * pi/180) ).'; % "b" formation
        B = [epsilon+real(b) epsilon-real(b) epsilon+imag(b) epsilon-imag(b)].';
        
        if( ( A*x - B ) >= 0 )
            G = [ G Theta_Grid(i)];
        end 
    end
end

%% Finda E-Matrix
function [E] = findE(ka,gamma, G,x)
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