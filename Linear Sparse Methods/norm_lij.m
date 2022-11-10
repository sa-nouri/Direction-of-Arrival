function f = norm_lij(A, i, j)
 
    for k = 1:length(A)
        f(k) = norms(A(k,:), i);
    end
   % x = norms(f,j);
end
