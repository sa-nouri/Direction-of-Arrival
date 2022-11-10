function x =  l1_svd(y,A_g ,Num_Source ,Num_Sample)
    fprintf(1, 'l2,1 norm optimization for multiple snapshots\n');

    Y = sv(y, Num_Source,Num_Sample);
    lambda = 0.5;
    cvx_begin
        variable S_sv(length(A_g), Num_Source)
        minimize ( power(2,norm((Y - A_g * S_sv ) , 'fro')) + lambda * sum(norms(S_sv)));
    cvx_end
    x = S_sv;
end