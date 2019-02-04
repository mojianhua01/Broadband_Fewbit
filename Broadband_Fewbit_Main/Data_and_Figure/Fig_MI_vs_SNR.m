% plot simulation results
% y: Rate
% x: SNR in dB
% The figure looks bad!

close all;
%load results_11032015.mat
% load results_QPSK_training_11032015.mat
% load results_QPSK_Nb_2048_07082016.mat
load results_Chu_Np_2048_08232016.mat
load results_VAMP_Chu_Np_2048_01012017.mat
% load results_QPSK_07082016.mat
% load results_Chu_10db_Nb_2048_bit_1_inf_07112016.mat
Nd = 16;
Nt = 64;
Nr = 64;

Training_length_index = 1;
% SNR_index = 1;

% MSE_BG = zeros(length(SNRdB_vector), length(bit_vector), length(Nb_vector), Num_chan);

% bit_vector = (1:1:9);

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

% MSE versus SNR
figure,
plot1 = plot(SNRdB_vector, MI_BG_mean(:,:, Training_length_index), 'linewidth', 1,'color','blue');
set(plot1(1),'DisplayName','1-bit','Marker','o','LineStyle','--');
set(plot1(2),'DisplayName','2-bit','Marker','square','LineStyle','--');
set(plot1(3),'DisplayName','3-bit','Marker','diamond','LineStyle','--');
set(plot1(4),'DisplayName','4-bit','Marker','v','LineStyle','--');
set(plot1(5),'DisplayName','\infty-bit','Marker','^','LineStyle','--');

grid on;
hold on;
plot1 = plot(SNRdB_vector, MI_GM_mean(:,:, Training_length_index), 'linewidth', 1,'color','black');
set(plot1(1),'DisplayName','1-bit','Marker','o','LineStyle','-');
set(plot1(2),'DisplayName','2-bit','Marker','square','LineStyle','-');
set(plot1(3),'DisplayName','3-bit','Marker','diamond','LineStyle','-');
set(plot1(4),'DisplayName','4-bit','Marker','v','LineStyle','-');
set(plot1(5),'DisplayName','\infty-bit','Marker','^','LineStyle','-');

plot1 = plot(SNRdB_vector, MI_VAMP_mean(:,:, Training_length_index), 'linewidth', 1,'color','red');
set(plot1(1),'DisplayName','1-bit','Marker','o','LineStyle','-');
set(plot1(2),'DisplayName','2-bit','Marker','square','LineStyle','-');
set(plot1(3),'DisplayName','3-bit','Marker','diamond','LineStyle','-');
set(plot1(4),'DisplayName','4-bit','Marker','v','LineStyle','-');
set(plot1(5),'DisplayName','\infty-bit','Marker','^','LineStyle','-');

plot1=plot(SNRdB_vector, MI_pcsi_mean(:,:, Training_length_index), 'linewidth', 1,'color',...
'magenta','LineStyle','-.');
set(plot1(1),'DisplayName','1-bit','Marker','o','LineStyle','-.');
set(plot1(2),'DisplayName','2-bit','Marker','square','LineStyle','-.');
set(plot1(3),'DisplayName','3-bit','Marker','diamond','LineStyle','-.');
set(plot1(4),'DisplayName','4-bit','Marker','v','LineStyle','-.');
set(plot1(5),'DisplayName','\infty-bit','Marker','^','LineStyle','-.');

xlabel('SNR (dB)', 'fontsize',14)
ylabel('Mutual Information (bps/Hz)', 'fontsize',14)
legend('EM-BG-GAMP, 1-bit', 'EM-BG-GAMP, 2-bit','EM-BG-GAMP, 3-bit', 'EM-BG-GAMP, 4-bit','EM-BG-GAMP, \infty-bit', ...
    'EM-GM-GAMP, 1-bit','EM-GM-GAMP, 2-bit', 'EM-GM-GAMP, 3-bit','EM-GM-GAMP, 4-bit', 'EM-GM-GAMP, \infty-bit', ...
    'EM-GM-VAMP, 1-bit','EM-GM-VAMP, 2-bit', 'EM-GM-VAMP, 3-bit','EM-GM-VAMP, 4-bit', 'EM-GM-VAMP, \infty-bit', ...
    'Perfect CSI, 1-bit', 'Perfect CSI, 2-bit','Perfect CSI, 3-bit', 'Perfect CSI, 4-bit','Perfect CSI, \infty-bit',...
    'Location','best')


%% plot 4-bit and infinite-bit results in a single figure
figure, 
plot_LS = plot(SNRdB_vector, MI_LS_mean(:,[4 5], Training_length_index), 'linewidth', 1,'Marker','d','color','[0 0.447 0.741]');
hold on;
grid on;
plot(SNRdB_vector, MI_LMMSE_mean(:,[4 5], Training_length_index), 'linewidth', 1,'Marker','v','color','[0.85 0.325 0.098]');
plot(SNRdB_vector, MI_BPDN_mean(:,[4 5], Training_length_index), 'linewidth', 1,'Marker','^','color','[0.929 0.694 0.125]');
plot(SNRdB_vector, MI_BG_mean(:,[4 5], Training_length_index), 'linewidth', 1,'Marker','o','color','[0.494 0.184 0.556]');
plot(SNRdB_vector, MI_GM_mean(:,[4 5], Training_length_index), 'linewidth', 1,'Marker','s','color','[0.466 0.674 0.188]');
xlabel('SNR (dB)', 'fontsize',14)
ylabel('Mutual Information (bps/Hz)', 'fontsize',14)
legend('LS, 1-bit', 'LS, \infty-bit','ALMMSE, 1-bit', 'ALMMSE, \infty-bit','BPDN, 1-bit', 'BPDN, \infty-bit',...
    'BG, 1-bit', 'BG, \infty-bit','GM, 1-bit', 'GM, \infty-bit', 'Location','best')