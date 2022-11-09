%% Capon

function [F_Capon] = Capon(R, d)
    angle = 0:0.01:180;
    a = ( 1/sqrt(length(d)) ) * exp( -1i * d' * pi .* cosd(angle));
    F_Capon = 1./sum(conj(a) .* (inv(R) * a));
    %F_Capon = diag(1./(a' * inv(R) * a));
end

%%
