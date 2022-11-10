function B = sracv(R , Num_Sample, d) 

    Theta_g = 0:0.1:180 ;
    A_g = exp( -1i .* d' .* pi .* cosd(Theta_g));


    fprintf(1, 'L1-norm --- Sparse Reconstruction of Array Covariance Vector (SRACV) \n');

    sigma = 1 ; M = length(d) ;
    Rn = (sigma ^2 ) * eye(M) ; 
    W = ones(length(Theta_g), 1);
    p = 0.0001;
    Etha = chi2inv(1-p, M^2);
    y = vec(sqrt(Num_Sample) * R^(-1/2) * (R-Rn) * R^(-1/2));
    Phi = sqrt(Num_Sample) * kron((R^(-1/2)).' ,(R^(-1/2) * A_g));

    cvx_begin  
        variable B(length(Theta_g),M) complex
        variable g
        cvx_solver sedumi
        minimize (g) 
        subject to 
            W' * (norms(B,2,2)) <= g ;
            norms(y - (Phi * vec(B))) <= sqrt(Etha);
    cvx_end


end
