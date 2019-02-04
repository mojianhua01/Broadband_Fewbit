% function [H_scaled, H_v]= func_mmWave_Channel_ULA(BS_ant, MS_ant,Num_cluster, Num_ray, Angle_spread)
%Generate

% %% setup for test
clear;
clc;
close all;
% for channel_index=1:1:100
BS_ant = 64;
MS_ant = 64;
Num_cluster = 2;
Num_ray = 4;
Angle_spread = pi/24; % pi/24;

% Simulation parameters
BS_ant_index = 0:1:BS_ant - 1;
MS_ant_index = 0:1:MS_ant - 1;
H = zeros(MS_ant,BS_ant);
KD = pi;  % Assuming: K=2pi/lambda, D=lambda/2

% Channel generation
for L = 1:1:Num_cluster
    alpha(1,L) = 1;%(sqrt(1/2)*(randn(1,1)+1j*randn(1,1))); % gain
    AoA_phi_mean = pi*rand(1,1) - pi/2;   % azimuthal angle
    AoD_phi_mean = pi*rand(1,1) - pi/2;
    
    for ii = 1:1:Num_ray
%         alpha( ii + (L-1) * Num_ray) = randn(1,2)*[1;1i]/sqrt(2); % gain
        alpha( ii + (L-1) * Num_ray) = 1;
        AoA_offset = (-1 + 2* rand(1,1))*Angle_spread;
        AoD_offset = (-1 + 2* rand(1,1))*Angle_spread;
        
        AoA_phi = AoA_phi_mean + AoA_offset;
        AoD_phi = AoD_phi_mean + AoD_offset;
        
        BS_steering_temp = sqrt(1/BS_ant)*exp(1j * KD * BS_ant_index * sin(AoD_phi));
        
        BS_steering(:,ii+(L-1)*Num_ray) = BS_steering_temp(:);
        
        MS_steering_temp = sqrt(1/MS_ant)*exp(1j * KD * MS_ant_index * sin(AoA_phi));
        
        MS_steering(:,ii+(L-1)*Num_ray) = MS_steering_temp(:);
        
        H = H + MS_steering * BS_steering' * exp(-AoA_offset) * exp(-AoD_offset) * alpha(L);
    end;
end

H_scaled = H.*sqrt( 1/sum(sum(abs(H).^2))) * sqrt(BS_ant*MS_ant);
% H_scaled = H;

H_v = 1/sqrt(BS_ant * MS_ant) * dftmtx(MS_ant)'*H_scaled*dftmtx(BS_ant);

H_v = H_v./norm(H_v, 'fro')*sqrt(BS_ant * MS_ant);
    norm(H_v, 'fro')

% Hv_all(:,:,channel_index) = H_v;
% end;
% save H_64_16_2paths.mat Hv_all;
%
figure(1)
subplot(1,2,1)
bar3(abs(H_scaled))
xlabel('TX antenna j')
ylabel('RX antenna i')
set(gca, 'FontSize',16)
title('Spatial Domain', 'FontSize', 20)
zlabel('$|\mathbf{H}_{i,j}|$', 'Interpreter', 'latex');

subplot(1,2,2)
bar3(abs(H_v))
% stem3(1:1:BS_ant ,1:1:MS_ant, (abs(H_v)))
set(gca, 'FontSize',16)
grid on;
xlabel('TX angle j')
ylabel('RX agnle i')
title('Angular Domain', 'FontSize', 20)
zlabel('$|\mathbf{X}_{i,j}|$', 'Interpreter', 'latex');

figure(2)
bar3(abs(H_scaled))
xlabel('TX antenna j')
ylabel('RX antenna i')
set(gca, 'FontSize',14)
zlabel('$|\mathbf{H}_{i,j}|$', 'Interpreter', 'latex');

figure(3)
bar3(abs(H_v))
xlabel('TX beam j')
ylabel('RX beam i')
set(gca, 'FontSize',14)
zlabel('$|\mathbf{X}_{i,j}|$', 'Interpreter', 'latex');

figure,
contour(1:1:BS_ant ,1:1:MS_ant, log10(abs(H_v)))
set(gca, 'FontSize',16)
grid on;
xlabel('TX beam #')
ylabel('RX beam #')


y = sort(abs(H_v(:)), 'descend');
figure,
loglog(y);
% xlim([0 300])
ylim([0 10 * y(1)])
temp = norm(alpha) *( 1 + log(BS_ant * MS_ant)); 
temp = sqrt(6)/pi * norm(y);
hold on;  plot( temp * ( (1:1:BS_ant*MS_ant)).^(-1), 'r--')
% plot(  y(1)* exp(- lambertw((1:1:BS_ant*MS_ant)/4)), 'k--')
ylabel('sorted |x|')
xlabel('index')
set(gca, 'fontsize', 12)
title(['#clusters=', num2str( Num_cluster), ', #subpaths=', num2str(Num_ray) ] )

% figure, cdfplot(y)
% hold on, plot(1- 1./(1:1:20).*log(1./(1:1:20)))