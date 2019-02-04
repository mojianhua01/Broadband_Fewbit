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

index_channel = 1;
for index_delay = 1:1:Nd
    
    % virtual channel model or angular domain channel model
    %         Gtrue = func_mmWave_Channel_UPA(Nt_W, Nt_H, Nr_W, Nr_H, ...
    %             Angle_spread_W, Angle_spread_H, Num_cluster, Num_ray, Narrowband_flag, Nd);
    Gtrue = Hv_all(:,:,index_delay,index_channel);
    Gtrue_norm(1, index_delay) = norm(Gtrue(:))^2;
end;

figure,
stem(Gtrue_norm)
xlabel('Delay tap $\ell$', 'fontsize', 14, 'interpreter', 'latex')
ylabel('$$\|\mathbf{H}[\ell]\|_F^2$$','fontsize',14,'interpreter','latex')