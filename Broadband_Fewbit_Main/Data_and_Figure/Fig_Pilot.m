% plot simulation results
% y: Rate
% x: SNR in dB
% The figure looks bad!

close all;
%load results_11032015.mat
% load results_QPSK_training_11032015.mat
% load results_QPSK_Nb_2048_07082016.mat

% load results_QPSK_07082016.mat
% load results_Chu_10db_Nb_2048_bit_1_inf_07112016.mat
Nd = 16;
Nt = 64;
Nr = 64;

Training_length_index = 1;

%% GM
% MSE versus SNR
figure,
% load results_GM_Gaussian_Np_2048_07122016.mat
load results_VAMP_Gaussian_Np_2048_01172017.mat
plot1 = plot(SNRdB_vector, NMSE_VAMP_dB(:,[1 2 3 4 5], Training_length_index), 'linewidth', 1, 'Marker', '^');
set(plot1(1),'color', '[0 0.447 0.741]');
set(plot1(2),'color', '[0.85 0.325 0.098]');
set(plot1(3),'color', '[0.929 0.694 0.125]');
set(plot1(4),'color', '[0.494 0.184 0.556]');
set(plot1(5),'color', '[0.466 0.674 0.188]');
xlabel('SNR (dB)', 'fontsize',14)
ylabel('Normalized MSE (dB)', 'fontsize',14)
grid on;
hold on;

% load results_GM_QPSK_Np_2048_07122016.mat
load results_VAMP_QPSK_Np_2048_01172017.mat
plot1 = plot(SNRdB_vector, NMSE_VAMP_dB(:,[1 2 3 4 5], Training_length_index), 'linewidth', 1, 'Marker', 'v');
set(plot1(1),'color', '[0 0.447 0.741]');
set(plot1(2),'color', '[0.85 0.325 0.098]');
set(plot1(3),'color', '[0.929 0.694 0.125]');
set(plot1(4),'color', '[0.494 0.184 0.556]');
set(plot1(5),'color', '[0.466 0.674 0.188]');

% load results_GM_Golay_Np_2048_09252016.mat
load results_VAMP_Golay_Np_2048_01172017.mat
plot1 = plot(SNRdB_vector, NMSE_VAMP_dB(:,[1 2 3 4 5], Training_length_index), 'linewidth', 1, 'Marker', 'd');
set(plot1(1),'color', '[0 0.447 0.741]');
set(plot1(2),'color', '[0.85 0.325 0.098]');
set(plot1(3),'color', '[0.929 0.694 0.125]');
set(plot1(4),'color', '[0.494 0.184 0.556]');
set(plot1(5),'color', '[0.466 0.674 0.188]');

% load results_GM_Chu_Np_2048_07122016.mat
load results_VAMP_Chu_Np_2048_01172017.mat
plot1 = plot(SNRdB_vector, NMSE_VAMP_dB(:,[1 2 3 4 5], Training_length_index), 'linewidth', 1, 'Marker', 's');
set(plot1(1),'color', '[0 0.447 0.741]');
set(plot1(2),'color', '[0.85 0.325 0.098]');
set(plot1(3),'color', '[0.929 0.694 0.125]');
set(plot1(4),'color', '[0.494 0.184 0.556]');
set(plot1(5),'color', '[0.466 0.674 0.188]');

% %% plot 1-bit and infinite-bit results in a single figure
% figure, 
% plot_LS = plot(SNRdB_vector, Rate_LS_mean(:,[1, 5], Training_length_index), 'linewidth', 1,'Marker','d','color','[0 0.447 0.741]');
% hold on;
% grid on;
% plot(SNRdB_vector, Rate_LMMSE_mean(:,[1, 5], Training_length_index), 'linewidth', 1,'Marker','v','color','[0.85 0.325 0.098]');
% plot(SNRdB_vector, Rate_BPDN_mean(:,[1, 5], Training_length_index), 'linewidth', 1,'Marker','^','color','[0.929 0.694 0.125]');
% plot(SNRdB_vector, Rate_BG_mean(:,[1, 5], Training_length_index), 'linewidth', 1,'Marker','o','color','[0.494 0.184 0.556]');
% plot(SNRdB_vector, Rate_GM_mean(:,[1, 5], Training_length_index), 'linewidth', 1,'Marker','s','color','[0.466 0.674 0.188]');
% xlabel('SNR (dB)', 'fontsize',14)
% ylabel('Normalized MSE (dB)', 'fontsize',14)
% legend('LS, 1-bit', 'LS, \infty-bit','LMMSE, 1-bit', 'LMMSE, \infty-bit','BPDN, 1-bit', 'BPDN, \infty-bit',...
%     'BG, 1-bit', 'BG, \infty-bit','GM, 1-bit', 'GM, \infty-bit', 'Location','best')