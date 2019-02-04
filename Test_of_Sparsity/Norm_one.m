clc;
clear;
% load H_UPA_64_64_16_4clusters.mat
% Nta = 8;  % azimuth
% Nte = 8;  % elevation
% Nt = Nta * Nte;
% Nra = 8;  % azimuth
% Nre = 8;  % elevation
% Nr = Nra * Nre;
% Nd = 16;  % Delay spread L

load H_UPA_16_4_16_2clusters.mat
Nta = 4;  % azimuth departure angle
Nte = 4;  % elevation departure angle
Nt = Nta * Nte;
Nra = 2;  % azimuth arrival angle
Nre = 2;  % elevation
Nr = Nra * Nre;
Nd = 16;  % Delay spread L

Num_chan = 100;
for index_channel = 1:1:Num_chan
    
    % virtual channel model or angular domain channel model
    %         Gtrue = func_mmWave_Channel_UPA(Nt_W, Nt_H, Nr_W, Nr_H, ...
    %             Angle_spread_W, Angle_spread_H, Num_cluster, Num_ray, Narrowband_flag, Nd);
    Gtrue = Hv_all(:,:,:,index_channel);
%     Gtrue = Hv_all(:,:,:,3);
    %     Gtrue = Gtrue/norm(Gtrue(:));
    gtrue = reshape(Gtrue, [],1);
    
    Gnorm_1(index_channel) = norm(gtrue,1);
    
    Gnorm_2(index_channel) = norm(gtrue,2);
end;

subplot(2,1,1)
plot(Gnorm_1);
hold on;
plot(Gnorm_2 * log(Nt*Nr*Nd));
legend('L1-norm','L2-norm * log N')
xlabel('samples')
mean(Gnorm_1./Gnorm_2/log(Nt*Nr*Nd))
title('Channel with 4 clusters')
subplot(2,1,2);
plot(Gnorm_1./Gnorm_2/log(Nt*Nr*Nd))
h = legend('$\frac{L1-norm}{L2-norm * \log N}$');
set(h,'Interpreter','latex')

sqrt(2*6/pi^2)