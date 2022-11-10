function Y = sv(y, Num_Source, Num_Sample) 
    
    [U, L, V] = svd(y);
    D = [eye(Num_Source) zeros(Num_Source, (Num_Sample - Num_Source))];
    Y = U * L * D';
    %Y =  y * V * D' ; 
end