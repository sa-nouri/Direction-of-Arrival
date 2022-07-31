%%
clc ; clear ; close all

%% Angle Estimation 
ItrMax = 100; Num_Source = 4 ; Num_Sample = 1e3;  SNR_dB = [10 10 10 10 ] ; M = 10 ;  % Num of Antennas 
Tetha = [ 20 ; 40 ; 60 ; 85 ] ; %rand(Num_Source,1) .* 180 / pi ;

error_beam = zeros(ItrMax,Num_Source);
error_capon = zeros(ItrMax,Num_Source);

for i = 1:ItrMax

    [R,A] = Correlation(M,Tetha,Num_Source,Num_Sample,SNR_dB);
    f_beam = Beamformer(R,M);
    f_capon = Capon(R,M);
    
 % finding peaks

    [Pks_beam,loc_beam ] = findpeaks(abs(f_beam),'NPeaks',Num_Source, 'SortStr','descend');
    [Pks_capon,loc_capon ] = findpeaks(abs(f_capon),'NPeaks',Num_Source, 'SortStr','descend');
    angle = 0:0.01:180;
    angle_beam = angle(loc_beam) ;  
    angle_capon = angle(loc_capon) ;

%% Angel Error Calculation
    for j = 1 : Num_Source
        error_beam(i,j) = Tetha(j,1) - angle_beam(1,j) ;
        error_capon(i,j) = Tetha(j,1) - angle_capon(1,j) ; 
    end
    
end


%% RMSE

rmse_beam = sqrt(sum((error_beam.^2))); 
rmse_capon = sqrt(sum(error_capon.^2));

% plot((rmse_beam),'b');
% grid on ;
% hold on ;
% plot((rmse_capon),'r');
% title(' RMSE  ' ) ;

%%




