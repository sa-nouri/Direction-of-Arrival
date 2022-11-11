%% ------------  Dual Optimization for Davis Transform  ******

close all; clear; clc;

%% Covariance Calculation & Initializing

Num_Sample = 1; SNR_dB = 100;   SNR = 10 .^ (SNR_dB/10) ;
M = 11; d = 0: M-1;  Theta = [0];  Phi = [24]';
Num_Source = length(Theta); gamma = 2*pi/M * d' * 180 /pi;

ka = pi/2;
[R, ~] = CirCov(d, Theta, Phi, Num_Source, Num_Sample, SNR_dB, ka);

%% Transformation UCA in Signal Environment

% Acording to algorithm which is attached in the paper
% Coding is based on Steps 

%% Step (1) --- Initialization

N = M; 
%x = rand(2*N,1) * 0;  % use inital value of paper
delta = 1e-5; I = 8; 

R = eye(N)/2; 
Q = [real(R) -imag(R) ; imag(R) real(R)];
Epsilon = [];
epsilon = 0.01;

%%

kk0 = 1: (2^4) + 1;
Theta_kk = (360)/(length(kk0) - 1) .* (kk0 - 1) ;
lambda_0 = rand(4, 1);

temp = zeros( 2 * N, 1);
for j = 1 : length(kk0)
    a = exp(1i .* ka .* cosd((Theta_kk(j) - gamma))).';
    A = [real(a) -imag(a) ; -real(a) imag(a) ; imag(a) real(a); - imag(a) -real(a)] ;
    temp = temp + (A.' * lambda_0);
end

x = - inv(Q) * temp;
%% find M0

M0_temp = (0: (N - 1)/2);
M0 = [];
for i = 1: length(M0_temp)
    temp = besselj(M0_temp - N, ka)/besselj(M0_temp, ka);
    if ( (temp <delta) && (M0_temp(i) <= (M-1))/2 )
        M0 = [M0 M0_temp];
    end
end
M0 = max(M0);

%% Itertation
i = 0;
temp = delta;

XX = [];

for l = -M0 :M0
    while((i < I) || (temp > delta))
        i = i + 1;
        K = (1: (2^(4+i)) + 1) ;
        Theta_Grid = (360)/(length(K)-1) .* ( K - 1) ;
        [Theta_ach] = findG(ka, gamma, x, Theta_Grid, Epsilon, epsilon, l); % G is not Create.
        if ~isempty(Theta_ach)
            [X, lambda] = DQSIP8(ka, gamma, Theta_ach, Q, epsilon, l); % (P8k) Solution
            temp = abs((X.' * Q * X) - ( x.' * Q * x)) ;
            [E] = findE(ka, gamma, Theta_ach, X, epsilon, l);
            Epsilon = E;
            x = X;
        end
    end

    [X, lambda, Theta_Opt] = DualQSIP7(ka, gamma, K, N, Q, l); % (P7k) Solution
    XX = [XX ; X];     
end

%% Circular Covariance Matrix
function [R, A] = CirCov(d, Tetha, Phi, Num_Source, Num_Sample, SNR_dB, ka)
    
    M = length(d);
    gamma = 2 * pi/M * d' * 180 /pi;
    A = exp(1i .* ka .* sind((Phi' - gamma)).*cosd(Tetha)  );
    SNR = 10 .^ (SNR_dB/10) ; 
    
    noise  = (randn(length(d), Num_Sample) + ( 1i * randn(length(d), Num_Sample))) * (sqrt(1/2));
    signal = (randn(Num_Source, Num_Sample) + ( 1i * randn(Num_Source, Num_Sample))) .* (sqrt(SNR(:,1)/2));
    
    x = A * signal + noise;
    R = 1/M * (x * x');
    R = R/length(signal);
end

%% Problem (6) Solution ----- Quadratic Semi-Infinite Programming
function [x] = QSIP6(ka, gamma, N, Theta_Grid, Q)
    
    epsilon = 1e-5;
    for i = 1: length(Theta_Grid)
        a = exp(1i .* ka .* cosd((Theta_Grid(i) - gamma))).'; % gamma value
        A = [real(a) -imag(a) ; -real(a) imag(a) ; imag(a) real(a); - imag(a) -real(a)];
        b = exp(1i .* ( Theta_Grid(i) * pi/180)).'; % "b" formation
        B = [epsilon+real(b) epsilon-real(b) epsilon+imag(b) epsilon-imag(b)].';
        cvx_begin
            variable x(2*N, 1)
            minimize(1/2 * x.' * Q * x)
            subject to
                (A * x - B ) <= zeros(length(B), 1);
        cvx_end
    end 

end

%% Problem (7) Solution ---- Dual Parametirization ---- Quadratic Semi-Infinite Programming
function [x, lambda, Theta_Opt] = DualQSIP7(ka, gamma, K, N, Q, M0)
    cvx_begin
        variable x(2 * N, 1)
        variable lambda(4, length(K))
        variable Theta_Opt(length(K), 1)
        b = exp(1i .* M0 *( Theta_Opt * pi/180) );
        B = [epsilon + real(b); epsilon - real(b); epsilon + imag(b); epsilon - imag(b)].';            
        minimize((1/2 * x.' * Q * x) + trace(B * lambda))
        subject to 
                for j = 1 :length(K)
                    a = exp(1i .* ka .* cosd((Theta_Opt(j) - gamma))).';
                    A = [real(a) -imag(a) ; -real(a) imag(a) ; imag(a) real(a); - imag(a) -real(a)] ;
                    temp = temp + (A.' * lambda(:,j)) ;
                end
            (Q * x) + temp == zeros(length(Q*x), 1);
    cvx_end

end

%% Problem (8) Solution ---- Dual Quadratic Semi-Infinite Programming
function [x, lambda] = DQSIP8(ka, gamma, Theta_Grid, Q, epsilon, M0)
    
    k = length(Theta_Grid) ;
    temp = 0;
    cvx_begin
        variable x(2*(length(gamma)), 1)
        variable lambda(4, k)

        b = exp(1i .* M0 *( Theta_Grid * pi/180));
        B = [epsilon+real(b); epsilon- real(b); epsilon+imag(b); epsilon-imag(b)].';            
        minimize((1/2 * x.' * Q * x) + trace(B * lambda)); 
        
        subject to
            for j = 1 : k
                a = exp(1i .* ka .* cosd((Theta_Grid(j) - gamma))).';
                A = [real(a) -imag(a) ; -real(a) imag(a) ; imag(a) real(a); - imag(a) -real(a)].' ;
                temp = temp + sum(A * lambda) ; % sum should be trace but there is error % Size of A is not possible to have Trace
            end
            (Q * x) + temp == zeros(length(Q*x), 1);
    
    cvx_end

end

%% Find G-Matrix 
function [G] = findG(ka, gamma, x, Theta_Grid, Epsilon, epsilon, M0)
    
    G = Epsilon;
    for i = 1: length(Theta_Grid)
        a = exp(1i * ka * cosd((Theta_Grid(i) - gamma))).'; % gamma value & 
        A = [real(a) -imag(a) ; -real(a) imag(a) ; imag(a) real(a); -imag(a) -real(a)] ;
        b = exp( 1i .* M0 *( Theta_Grid(i) * pi/180)); % "b" formation
        B = [epsilon+real(b) epsilon-real(b) epsilon+imag(b) epsilon-imag(b)].';
        
        if (all(( A*x - B ) >= 0)) %% There is Problem : 2 elements of matrix are Positve and others are Negative. :))
            G = [G Theta_Grid(i)];
        end 
    end
    G = unique(G);
end

%% Finda E-Matrix
function [E] = findE(ka, gamma, G, x, epsilon, M0)
    E = [];
    for i = 1: length(G)
        a = exp(1i .* ka .* cosd((G(i) - gamma))).'; % gamma value & 
        A = [real(a) -imag(a) ; -real(a) imag(a) ; imag(a) real(a); - imag(a) -real(a)] ;
        b = exp( 1i .* M0 *  ( G(i) * pi/180) ); % "b" formation
        B = [epsilon+real(b) epsilon-real(b) epsilon+imag(b) epsilon-imag(b)].';
        
        if((A*x - B) == 0)
            E = [ E G(i)];
        end
    end
end