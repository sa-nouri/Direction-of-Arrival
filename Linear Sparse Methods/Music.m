%% MUSIC

function [F_Music] = Music(R,d,Num_Source)
    [Betha , lambda ] = eig(R);
    [~ , lambda_index ] = sort(diag(lambda),'descend');
    Betha = Betha(:,lambda_index); 
    U = Betha(:,Num_Source+1:end);
    angle = -90:0.1:90 ;
    a = ( 1/sqrt(length(d)) ) * exp( -1i * d' * pi .* sind(angle));
    F_Music =(1./sum(conj(a) .* (U * U' * a)));
    
end

%%