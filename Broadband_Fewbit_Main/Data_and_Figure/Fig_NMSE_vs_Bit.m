% plot simulation results
% y: NMSE
% x: Training Length

close all;
%load results_11032015.mat
% load results_QPSK_training_11032015.mat
% load results_QPSK_Nb_2048_07082016.mat
% load results_Chu_Nb_2048_07082016.mat
% load results_QPSK_07082016.mat
% load results_Chu_10db_Np_2048_bit_1_inf_08162016.mat
% load results_VAMP_Chu_10db_Np_2048_bit_1_inf_01012017.mat

load results_Chu_10db_Np_2048_bit_1_inf_02252017.mat
Nd = 16;
Nt = 64;
Nr = 64;

Training_length_index = 1;
SNR_index = 1;

% MSE_BG = zeros(length(SNRdB_vector), length(bit_vector), length(Nb_vector), Num_chan);

bit_vector = (1:1:9);
NMSE_BG_VAMP_dB = Results.NMSE.NMSE_BG_VAMP_dB;
NMSE_GM_VAMP_dB = Results.NMSE.NMSE_GM_VAMP_dB;
NMSE_BG_GAMP_dB = Results.NMSE.NMSE_BG_GAMP_dB;
NMSE_GM_GAMP_dB = Results.NMSE.NMSE_GM_GAMP_dB;

NMSE_LS_dB = Results.NMSE.NMSE_LS_dB;
NMSE_LMMSE_dB = Results.NMSE.NMSE_LMMSE_dB;
NMSE_BPDN_dB = Results.NMSE.NMSE_BPDN_dB;
NMSE_QIHT_dB = Results.NMSE.NMSE_QIHT_dB;

% % MSE versus bit
figure,
plot_LS = plot(bit_vector, NMSE_LS_dB(SNR_index, :, Training_length_index), 'linewidth', 1,'Marker','d','color','[0 0.447 0.741]');
hold on;
grid on;
plot(bit_vector, NMSE_LMMSE_dB(SNR_index, :, Training_length_index), 'linewidth', 1,'Marker','x','color','[0.85 0.325 0.098]');
plot(bit_vector, NMSE_BPDN_dB(SNR_index, :, Training_length_index), 'linewidth', 1,'Marker','v','color','[0.929 0.694 0.125]');
plot(bit_vector, NMSE_QIHT_dB(SNR_index, :, Training_length_index), 'linewidth', 1,'Marker','<','color','[0.494 0.184 0.556]');
plot(bit_vector, NMSE_BG_GAMP_dB(SNR_index, :, Training_length_index), 'linewidth', 1,'Marker','o','color','[0.301 0.745 0.933]');
plot(bit_vector, NMSE_GM_GAMP_dB(SNR_index, :, Training_length_index), 'linewidth', 1,'Marker','s','color','[0.635 0.078 0.184]');
plot(bit_vector, NMSE_BG_VAMP_dB(SNR_index, :, Training_length_index), 'linewidth', 1,'Marker','^','color','[0.466 0.674 0.188]');
plot(bit_vector, NMSE_GM_VAMP_dB(SNR_index, :, Training_length_index), 'linewidth', 1,'Marker','+');

xlabel('ADC Resolution (bit)', 'fontsize',14)
ylabel('Normalized MSE (dB)', 'fontsize',14)
legend('LS', 'ALMMSE','SPGL1','QIHT','EM-BG-GAMP','EM-GM-GAMP','EM-BG-VAMP','EM-GM-VAMP','Location','best')
xticks([1 2 3 4 5 6 7 8 9])
xticklabels({'1','2','3','4','5','6','7','8', '\infty'})
