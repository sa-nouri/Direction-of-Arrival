function [x] = BPDN(y, A)

    res = y ; % residual( y - Ax ) 
    norm_y = norm(y) ;
    y_temp = zeros(length(y), 1);

    %%shrink = @(z) sign(z) .* max(0,abs(z)-1);
    alpha = 1;

    max_it = 3e3; 
    max_time = 1e3;
    stepsize = 1e-3;    %1/max(abs(y.'*A));
    error = 1e-4;
    sigma = 1;

    for i = 1:max_it
        y_temp = y_temp + stepsize * res; % for k=1, (stepsize + 1/max(abs(b.'*A))) can be used as stepsize ---- Y update
     
        x = A' * y_temp ;
        %  x = alpha *  shrink(y .'* A ).';
        Ax = A * x;
        res = y - Ax; % res will be used in next y-update ----- X update
        diff = norm(res) + norm(x, 1) ----- Lasso 
        % diff = sum(norm(x, 2) .* W )  ----- Yall-1 
        % diff =  norm(x) && norm(x, 1) ---  MMV
        % diff = norm(res, 2) ---- CVX
         
        if norm(diff) < ( error); break; end 
    end
end

