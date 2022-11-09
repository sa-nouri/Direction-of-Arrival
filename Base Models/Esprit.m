% Calcluate Direction of Arrival with ESPRIT

function Esprit_Doas = Esprit(R, d, Num_Source)
    
    [Betha,lambda] = eig(R);
    [~,lambda_index] = sort(diag(lambda),1,'descend');
    Betha = Betha(:,lambda_index);
    
    U_signal = Betha(:,1:Num_Source);
    U_noise = Betha (:,Num_Source+1:length(d));
    
    Phi = linsolve(U_signal(1:length(d)-1,:),U_signal(2:length(d),:));
    Esprit_Doas = 90 - (asind(-angle(eig(Phi))/(pi)));

end
