% plot simulation results
% y: Rate
% x: bit

close all;
%load results_11032015.mat
% load results_QPSK_training_11032015.mat
% load results_QPSK_Nb_2048_07082016.mat
% load results_Chu_Nb_2048_07082016.mat
% load results_QPSK_07082016.mat
% load results_Chu_10db_Np_2048_bit_1_inf_08232016.mat
% load results_VAMP_Chu_10db_Np_2048_bit_1_inf_01012017.mat

load results_Chu_10db_Np_2048_bit_1_inf_02252017.mat
Nd = 16;
Nt = 64;
Nr = 64;

MI_BG_VAMP_mean = Results.MI.MI_BG_VAMP_mean;
MI_GM_VAMP_mean = Results.MI.MI_GM_VAMP_mean;
MI_BG_GAMP_mean = Results.MI.MI_BG_GAMP_mean;
MI_GM_GAMP_mean = Results.MI.MI_GM_GAMP_mean;

MI_LS_mean = Results.MI.MI_LS_mean;
MI_LMMSE_mean = Results.MI.MI_LMMSE_mean;
MI_BPDN_mean = Results.MI.MI_BPDN_mean;
MI_QIHT_mean = Results.MI.MI_QIHT_mean;

MI_pcsi_mean = Results.MI.MI_pcsi_mean;

Training_length_index = 1;
SNR_index = 1;

% MSE_BG = zeros(length(SNRdB_vector), length(bit_vector), length(Nb_vector), Num_chan);

bit_vector = (1:1:9);

% Nc = 10000;
% Np = 2048;
% prefactor = (Nc-Np)/Nc; 
% Rate_LS_mean = prefactor * Rate_LS_mean;
% Rate_LMMSE_mean = prefactor * Rate_LMMSE_mean;
% Rate_BPDN_mean = prefactor * Rate_BPDN_mean;
% Rate_BG_mean = prefactor * Rate_BG_mean;
% Rate_GM_mean = prefactor * Rate_GM_mean;
% Rate_QIST_mean = prefactor * Rate_QIST_mean;
% Rate_QIHT_mean = prefactor * Rate_QIHT_mean;

% % MSE versus bit
figure,
plot_LS = plot(bit_vector, MI_LS_mean(SNR_index, :, Training_length_index), 'linewidth', 1,'Marker','d','color','[0 0.447 0.741]');
hold on;
grid on;
plot(bit_vector, MI_LMMSE_mean(SNR_index, :, Training_length_index), 'linewidth', 1,'Marker','x','color','[0.85 0.325 0.098]');
plot(bit_vector, MI_BPDN_mean(SNR_index, :, Training_length_index), 'linewidth', 1,'Marker','v','color','[0.929 0.694 0.125]');
plot(bit_vector, MI_QIHT_mean(SNR_index, :, Training_length_index), 'linewidth', 1,'Marker','<','color','[0.494 0.184 0.556]');
plot(bit_vector, MI_BG_GAMP_mean(SNR_index, :, Training_length_index), 'linewidth', 1,'Marker','o','color','[0.301 0.745 0.933]');
plot(bit_vector, MI_GM_GAMP_mean(SNR_index, :, Training_length_index), 'linewidth', 1,'Marker','s','color','[0.635 0.078 0.184]');
plot(bit_vector, MI_BG_VAMP_mean(SNR_index, :, Training_length_index), 'linewidth', 1,'Marker','^','color','[0.466 0.674 0.188]');
plot(bit_vector, MI_GM_VAMP_mean(SNR_index, :, Training_length_index), 'linewidth', 1,'Marker','+');
plot(bit_vector, MI_pcsi_mean(SNR_index, :, Training_length_index), 'linewidth', 1,'Marker','p','color','[0 0.447 0.741]');

xlabel('ADC Resolution (bit)', 'fontsize',14)
ylabel('Mutual Information (bps/Hz)', 'fontsize',14)
legend('LS', 'ALMMSE','SPGL1','QIHT','EM-BG-GAMP','EM-GM-GAMP','EM-BG-VAMP','EM-GM-VAMP','Perfect CSI','Location','best')
xticks([1 2 3 4 5 6 7 8 9])
xticklabels({'1','2','3','4','5','6','7','8', '\infty'})
