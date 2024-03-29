%% Beamformer

function [F_Beam] = Beamformer(R, d)
    angle = 0:0.01:180;
    a = ( 1/sqrt(length(d)) ) * exp( -1i * d' * pi .* cosd(angle));
    F_Beam = sum(conj(a) .* (R * a));
    %F_Beam = diag( a' * R * a);
end

%%

