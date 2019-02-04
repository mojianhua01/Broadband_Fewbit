% plot simulation results
% y: Rate
% x: Training Length


close all;
% % load results_Chu_0dB_Nb_512_4096_07112016.mat
%
% load results_VAMP_Chu_10dB_Np_1024_5120_01012017.mat
% load results_GAMP_Chu_10dB_Np_1024_5120_01012017.mat
%
% Nd = 16;
% Nt = 64;
% Nr = 64;
%
% % Training_length_index = 1;
% SNR_index = 1;
%
% % Nc = 10000;
% % prefactor = (Nc - Np_vector)/Nc;
% % Rate_LS_mean = bsxfun(@times, prefactor',  Rate_LS_mean(SNR_index,:,:));
% % Rate_LMMSE_mean = prefactor * Rate_LMMSE_mean;
% % Rate_BPDN_mean = prefactor * Rate_BPDN_mean;
% % Rate_BG_mean = prefactor * Rate_BG_mean;
% % Rate_GM_mean = prefactor * Rate_GM_mean;
% % Rate_QIST_mean = prefactor * Rate_QIST_mean;
% % Rate_QIHT_mean = prefactor * Rate_QIHT_mean;
% %% BG
%
% % % MSE versus Np
% figure,
% plot1 = plot(Np_vector, ( reshape( Rate_BG_mean(SNR_index,:,:), length(bit_vector), []) )', 'linewidth', 1.5);
% set(plot1(1),'DisplayName','1-bit','Marker','o','LineStyle','-');
% set(plot1(2),'DisplayName','2-bit','Marker','square','LineStyle','-');
% set(plot1(3),'DisplayName','3-bit','Marker','diamond','LineStyle','-');
% set(plot1(4),'DisplayName','4-bit','Marker','v','LineStyle','-');
% % set(plot1(5),'DisplayName','\infty-bit','Marker','^','LineStyle','--');
% xlabel('Training length', 'fontsize',14)
% ylabel('Rate (bps/Hz)', 'fontsize',14)
% title('EM-BG-AMP algorithm')
% grid on;
% set(gca, 'XTick', Np_vector)
%
%
% %% GM
%
% % % Rate versus Np
% figure,
% plot1 = plot(Np_vector, ( reshape( Rate_GM_mean(SNR_index,:,:), length(bit_vector), []) )', 'linewidth', 1.5);
% set(plot1(1),'DisplayName','1-bit','Marker','o','LineStyle','-');
% set(plot1(2),'DisplayName','2-bit','Marker','square','LineStyle','-');
% set(plot1(3),'DisplayName','3-bit','Marker','diamond','LineStyle','-');
% set(plot1(4),'DisplayName','4-bit','Marker','v','LineStyle','-');
% % set(plot1(5),'DisplayName','\infty-bit','Marker','^','LineStyle','--');
% xlabel('Training length', 'fontsize',14)
% ylabel('Rate (bps/Hz)', 'fontsize',14)
% title('EM-GM-AMP algorithm')
% grid on;
% set(gca, 'XTick', Np_vector)
%
%
% % %% BPDN
% %
% % % % Rate versus Np
% % figure,
% % plot1 = plot(Np_vector, ( reshape( Rate_BPDN_mean(SNR_index,:,:), length(bit_vector), []) )', 'linewidth', 1.5);
% % set(plot1(1),'DisplayName','1-bit','Marker','o','LineStyle','-');
% % set(plot1(2),'DisplayName','2-bit','Marker','square','LineStyle','-');
% % set(plot1(3),'DisplayName','3-bit','Marker','diamond','LineStyle','-');
% % set(plot1(4),'DisplayName','4-bit','Marker','v','LineStyle','-');
% % % set(plot1(5),'DisplayName','\infty-bit','Marker','^','LineStyle','--');
% % xlabel('Training length', 'fontsize',14)
% % ylabel('Rate (bps/Hz)', 'fontsize',14)
% % title('BPDN algorithm')
% % grid on;
% % set(gca, 'XTick', Np_vector)
% %
% % %% LS
% %
% % % % Rate versus Np
% % figure,
% % plot1 = plot(Np_vector, ( reshape( Rate_LS_mean(SNR_index,:,:), length(bit_vector), []) )', 'linewidth', 1.5);
% % set(plot1(1),'DisplayName','1-bit','Marker','o','LineStyle','-');
% % set(plot1(2),'DisplayName','2-bit','Marker','square','LineStyle','-');
% % set(plot1(3),'DisplayName','3-bit','Marker','diamond','LineStyle','-');
% % set(plot1(4),'DisplayName','4-bit','Marker','v','LineStyle','-');
% % % set(plot1(5),'DisplayName','\infty-bit','Marker','^','LineStyle','--');
% % xlabel('Training length', 'fontsize',14)
% % ylabel('Rate (bps/Hz)', 'fontsize',14)
% % title('LS algorithm')
% % grid on;
% % set(gca, 'XTick', Np_vector)
% %
% %
% % %% LMMSE
% %
% % % % Rate versus Np
% % figure,
% % plot1 = plot(Np_vector, ( reshape( Rate_LMMSE_mean(SNR_index,:,:), length(bit_vector), []) )', 'linewidth', 1.5);
% % set(plot1(1),'DisplayName','1-bit','Marker','o','LineStyle','-');
% % set(plot1(2),'DisplayName','2-bit','Marker','square','LineStyle','-');
% % set(plot1(3),'DisplayName','3-bit','Marker','diamond','LineStyle','-');
% % set(plot1(4),'DisplayName','4-bit','Marker','v','LineStyle','-');
% % % set(plot1(5),'DisplayName','\infty-bit','Marker','^','LineStyle','--');
% % xlabel('Training length', 'fontsize',14)
% % ylabel('Rate (bps/Hz)', 'fontsize',14)
% % title('LMMSE algorithm')
% % grid on;
% % set(gca, 'XTick', Np_vector)
% %
% % %% QIHT
% %
% % % % Rate versus Np
% % figure,
% % plot1 = plot(Np_vector, ( reshape( Rate_QIHT_mean(SNR_index,:,:), length(bit_vector), []) )', 'linewidth', 1.5);
% % set(plot1(1),'DisplayName','1-bit','Marker','o','LineStyle','-');
% % set(plot1(2),'DisplayName','2-bit','Marker','square','LineStyle','-');
% % set(plot1(3),'DisplayName','3-bit','Marker','diamond','LineStyle','-');
% % set(plot1(4),'DisplayName','4-bit','Marker','v','LineStyle','-');
% % % set(plot1(5),'DisplayName','\infty-bit','Marker','^','LineStyle','--');
% % xlabel('Training length', 'fontsize',14)
% % ylabel('Rate (bps/Hz)', 'fontsize',14)
% % title('QIHT algorithm')
% % grid on;
% % set(gca, 'XTick', Np_vector)
%
% % %% QIST
% %
% % % % Rate versus Np
% % figure,
% % plot1 = plot(Np_vector, ( reshape( Rate_QIST_mean(SNR_index,:,:), length(bit_vector), []) )', 'linewidth', 1.5);
% % set(plot1(1),'DisplayName','1-bit','Marker','o','LineStyle','-');
% % set(plot1(2),'DisplayName','2-bit','Marker','square','LineStyle','-');
% % set(plot1(3),'DisplayName','3-bit','Marker','diamond','LineStyle','-');
% % set(plot1(4),'DisplayName','4-bit','Marker','v','LineStyle','-');
% % set(plot1(5),'DisplayName','\infty-bit','Marker','^','LineStyle','--');
% % xlabel('Training length', 'fontsize',14)
% % ylabel('Rate (bps/Hz)', 'fontsize',14)
% % title('QIST algorithm')
% % grid on;
% % set(gca, 'XTick', Np_vector)
%
%
% %% plot few-bit and infinite-bit results in a single figure
% figure,
% plot1 = plot(Np_vector, ( reshape( Rate_BG_mean(SNR_index,:,:), length(bit_vector), []) )', 'linewidth', 1, 'color', 'blue');
% set(plot1(1),'DisplayName','1-bit','Marker','o','LineStyle','--');
% set(plot1(2),'DisplayName','2-bit','Marker','square','LineStyle','--');
% set(plot1(3),'DisplayName','3-bit','Marker','diamond','LineStyle','--');
% set(plot1(4),'DisplayName','4-bit','Marker','v','LineStyle','--');
% set(plot1(5),'DisplayName','\infty-bit','Marker','^','LineStyle','--');
% hold on;
% grid on;
% plot1 = plot(Np_vector, ( reshape( Rate_GM_mean(SNR_index,:,:), length(bit_vector), []) )', 'linewidth', 1, 'color', 'black');
% set(plot1(1),'DisplayName','1-bit','Marker','o','LineStyle','-');
% set(plot1(2),'DisplayName','2-bit','Marker','square','LineStyle','-');
% set(plot1(3),'DisplayName','3-bit','Marker','diamond','LineStyle','-');
% set(plot1(4),'DisplayName','4-bit','Marker','v','LineStyle','-');
% set(plot1(5),'DisplayName','\infty-bit','Marker','^','LineStyle','-');
%
% plot1 = plot(Np_vector, ( reshape( Rate_VAMP_mean(SNR_index,:,:), 5 , []) )', 'linewidth', 1, 'color', 'red');
% set(plot1(1),'DisplayName','1-bit','Marker','o','LineStyle','-');
% set(plot1(2),'DisplayName','2-bit','Marker','square','LineStyle','-');
% set(plot1(3),'DisplayName','3-bit','Marker','diamond','LineStyle','-');
% set(plot1(4),'DisplayName','4-bit','Marker','v','LineStyle','-');
% set(plot1(5),'DisplayName','\infty-bit','Marker','^','LineStyle','-');
% % plot(Np_vector, ( reshape( MI_pcsi_mean(SNR_index,4,:), 1, []) )', 'linewidth', 1,'Marker','s','color','[0.635 0.078 0.184]');
%
% xlabel('Training Length', 'fontsize',14)
% ylabel('Achievable Rate (bps/Hz)', 'fontsize',14)
% legend('EM-BG-GAMP, 1-bit', 'EM-BG-GAMP, 2-bit','EM-BG-GAMP, 3-bit', 'EM-BG-GAMP, 4-bit','EM-BG-GAMP, \infty-bit', ...
%     'EM-GM-GAMP, 1-bit','EM-GM-GAMP, 2-bit', 'EM-GM-GAMP, 3-bit','EM-GM-GAMP, 4-bit', 'EM-GM-GAMP, \infty-bit', ...
%     'EM-GM-VAMP, 1-bit', 'EM-GM-VAMP, 2-bit','EM-GM-VAMP, 3-bit', 'EM-GM-VAMP, 4-bit','EM-GM-VAMP, \infty-bit',...
%     'Location','best')
% xlim([1000 7000])
% set(gca, 'XTick', Np_vector(1:2:end))
%
%
% load results_Chu_10dB_Np_1024_5120_09232016.mat
% %% plot 2-bit results in a single figure
% figure,
% bit_index = 2;
% plot_LS = plot(Np_vector, ( reshape( Rate_LS_mean(SNR_index,bit_index,:), 1, []) )', 'linewidth', 1,'Marker','d','color','[0 0.447 0.741]');
% hold on;
% grid on;
% plot(Np_vector, ( reshape( Rate_LMMSE_mean(SNR_index,bit_index,:), 1, []) )', 'linewidth', 1,'Marker','x','color','[0.85 0.325 0.098]');
% plot(Np_vector, ( reshape( Rate_BPDN_mean(SNR_index,bit_index,:), 1, []) )', 'linewidth', 1,'Marker','v','color','[0.929 0.694 0.125]');
% plot(Np_vector, ( reshape( Rate_QIHT_mean(SNR_index,bit_index,:), 1, []) )', 'linewidth', 1,'Marker','<','color','[0.494 0.184 0.556]');
%
% load results_GAMP_Chu_10dB_Np_1024_5120_01012017.mat
% plot(Np_vector, ( reshape( Rate_BG_mean(SNR_index,bit_index,:), 1, []) )', 'linewidth', 1,'Marker','o','color','[0.301 0.745 0.933]');
% plot(Np_vector, ( reshape( Rate_GM_mean(SNR_index,bit_index,:), 1, []) )', 'linewidth', 1,'Marker','s','color','[0.635 0.078 0.184]');
% plot(Np_vector, ( reshape( Rate_VAMP_mean(SNR_index,bit_index,:), 1, []) )', 'linewidth', 1,'Marker','^','color','[0.466 0.674 0.188]');
% %plot(Np_vector, ( reshape( MI_pcsi_mean(SNR_index,bit_index,:), 1, []) )', 'linewidth', 1,'Marker','+','color','[0 0.447 0.741]');
%
% xlabel('Training Length', 'fontsize',14)
% ylabel('Achievable Rate (bps/Hz)', 'fontsize',14)
% legend('LS, 2-bit','ALMMSE, 2-bit', 'SPGL1, 2-bit', 'QIHT, 2-bit', 'EM-BG-GAMP, 2-bit', ...
%     'EM-GM-GAMP, 2-bit', 'EM-GM-VAMP, 2-bit', 'Location','best')
% set(gca, 'XTick', Np_vector(1:2:end))


load results_Chu_10dB_VAMP_GAMP_Np_1024_5120_02252017.mat
SNR_index = 1;
Rate_BG_VAMP_mean = Results.Rate.Rate_BG_VAMP_mean;
Rate_GM_VAMP_mean = Results.Rate.Rate_GM_VAMP_mean;
Rate_BG_GAMP_mean = Results.Rate.Rate_BG_GAMP_mean;
Rate_GM_GAMP_mean = Results.Rate.Rate_GM_GAMP_mean;

NMSE_BG_VAMP_dB = Results.NMSE.NMSE_BG_VAMP_dB;
NMSE_GM_VAMP_dB = Results.NMSE.NMSE_GM_VAMP_dB;
NMSE_BG_GAMP_dB = Results.NMSE.NMSE_BG_GAMP_dB;
NMSE_GM_GAMP_dB = Results.NMSE.NMSE_GM_GAMP_dB;

Np_vector = Input.Np_vector;
bit_vector = Input.bit_vector;

%% plot few-bit and infinite-bit results in a single figure
figure,
plot1 = plot(Np_vector, ( reshape( Rate_BG_GAMP_mean(SNR_index,:,:), length(bit_vector), []) )', 'linewidth', 1, 'color', 'magenta');
set(plot1(1),'DisplayName','1-bit','Marker','o','LineStyle','--');
set(plot1(2),'DisplayName','2-bit','Marker','square','LineStyle','--');
set(plot1(3),'DisplayName','3-bit','Marker','diamond','LineStyle','--');
set(plot1(4),'DisplayName','4-bit','Marker','v','LineStyle','--');
set(plot1(5),'DisplayName','\infty-bit','Marker','^','LineStyle','--');
hold on;
grid on;
plot1 = plot(Np_vector, ( reshape( Rate_GM_GAMP_mean(SNR_index,:,:), length(bit_vector), []) )', 'linewidth', 1, 'color', 'red');
set(plot1(1),'DisplayName','1-bit','Marker','o','LineStyle','-');
set(plot1(2),'DisplayName','2-bit','Marker','square','LineStyle','-');
set(plot1(3),'DisplayName','3-bit','Marker','diamond','LineStyle','-');
set(plot1(4),'DisplayName','4-bit','Marker','v','LineStyle','-');
set(plot1(5),'DisplayName','\infty-bit','Marker','^','LineStyle','-');

plot1 = plot(Np_vector, ( reshape( Rate_BG_VAMP_mean(SNR_index,:,:), 5 , []) )', 'linewidth', 1, 'color', 'blue');
set(plot1(1),'DisplayName','1-bit','Marker','o','LineStyle','--');
set(plot1(2),'DisplayName','2-bit','Marker','square','LineStyle','--');
set(plot1(3),'DisplayName','3-bit','Marker','diamond','LineStyle','--');
set(plot1(4),'DisplayName','4-bit','Marker','v','LineStyle','--');
set(plot1(5),'DisplayName','\infty-bit','Marker','^','LineStyle','--');

plot1 = plot(Np_vector, ( reshape( Rate_GM_VAMP_mean(SNR_index,:,:), 5 , []) )', 'linewidth', 1, 'color', 'black');
set(plot1(1),'DisplayName','1-bit','Marker','o','LineStyle','-');
set(plot1(2),'DisplayName','2-bit','Marker','square','LineStyle','-');
set(plot1(3),'DisplayName','3-bit','Marker','diamond','LineStyle','-');
set(plot1(4),'DisplayName','4-bit','Marker','v','LineStyle','-');
set(plot1(5),'DisplayName','\infty-bit','Marker','^','LineStyle','-');

xlabel('Training Length', 'fontsize',14)
ylabel('Achievable Rate (bps/Hz)', 'fontsize',14)
legend({'EM-BG-GAMP, 1-bit', 'EM-BG-GAMP, 2-bit','EM-BG-GAMP, 3-bit', 'EM-BG-GAMP, 4-bit','EM-BG-GAMP, \infty-bit', ...
    'EM-GM-GAMP, 1-bit','EM-GM-GAMP, 2-bit', 'EM-GM-GAMP, 3-bit','EM-GM-GAMP, 4-bit', 'EM-GM-GAMP, \infty-bit', ...
    'EM-BG-VAMP, 1-bit', 'EM-BG-VAMP, 2-bit','EM-BG-VAMP, 3-bit', 'EM-BG-VAMP, 4-bit','EM-BG-VAMP, \infty-bit' ...
    'EM-GM-VAMP, 1-bit', 'EM-GM-VAMP, 2-bit','EM-GM-VAMP, 3-bit', 'EM-GM-VAMP, 4-bit','EM-GM-VAMP, \infty-bit'},...
    'Location','best')
xlim([1000 7000])
set(gca, 'XTick', Np_vector(1:2:end))


load results_Chu_Rate_Np_2bit_10dB_02252017.mat
Rate_BG_VAMP_mean = Results.Rate.Rate_BG_VAMP_mean;
Rate_GM_VAMP_mean = Results.Rate.Rate_GM_VAMP_mean;
Rate_BG_GAMP_mean = Results.Rate.Rate_BG_GAMP_mean;
Rate_GM_GAMP_mean = Results.Rate.Rate_GM_GAMP_mean;

Rate_LS_mean = Results.Rate.Rate_LS_mean;
Rate_LMMSE_mean = Results.Rate.Rate_LMMSE_mean;
Rate_BPDN_mean = Results.Rate.Rate_BPDN_mean;
Rate_QIHT_mean = Results.Rate.Rate_QIHT_mean;

%% plot 2-bit results in a single figure
figure,
bit_index = 1;
Rate_LS_mean(:,:,3) = 18;
Rate_LMMSE_mean(:,:,3) = 18;
plot_LS = plot(Np_vector, ( reshape( Rate_LS_mean(SNR_index,bit_index,:), 1, []) )', 'linewidth', 1,'Marker','d','color','[0 0.447 0.741]');
hold on;
grid on;
plot(Np_vector, ( reshape( Rate_LMMSE_mean(SNR_index,bit_index,:), 1, []) )', 'linewidth', 1,'Marker','x','color','[0.85 0.325 0.098]');
plot(Np_vector, ( reshape( Rate_BPDN_mean(SNR_index,bit_index,:), 1, []) )', 'linewidth', 1,'Marker','v','color','[0.929 0.694 0.125]');
plot(Np_vector, ( reshape( Rate_QIHT_mean(SNR_index,bit_index,:), 1, []) )', 'linewidth', 1,'Marker','<','color','[0.494 0.184 0.556]');

load results_GAMP_Chu_10dB_Np_1024_5120_01012017.mat
plot(Np_vector, ( reshape( Rate_BG_GAMP_mean(SNR_index,bit_index,:), 1, []) )', 'linewidth', 1,'Marker','o','color','[0.301 0.745 0.933]');
plot(Np_vector, ( reshape( Rate_GM_GAMP_mean(SNR_index,bit_index,:), 1, []) )', 'linewidth', 1,'Marker','s','color','[0.635 0.078 0.184]');
plot(Np_vector, ( reshape( Rate_BG_VAMP_mean(SNR_index,bit_index,:), 1, []) )', 'linewidth', 1,'Marker','^','color','[0.466 0.674 0.188]');
plot(Np_vector, ( reshape( Rate_GM_VAMP_mean(SNR_index,bit_index,:), 1, []) )', 'linewidth', 1,'Marker','+');
%plot(Np_vector, ( reshape( MI_pcsi_mean(SNR_index,bit_index,:), 1, []) )', 'linewidth', 1,'Marker','+','color','[0 0.447 0.741]');

xlabel('Training Length', 'fontsize',14)
ylabel('Achievable Rate (bps/Hz)', 'fontsize',14)
legend('LS, 2-bit','ALMMSE, 2-bit', 'SPGL1, 2-bit', 'QIHT, 2-bit', 'EM-BG-GAMP, 2-bit', ...
    'EM-GM-GAMP, 2-bit', 'EM-BG-VAMP, 2-bit','EM-GM-VAMP, 2-bit', 'Location','best')
set(gca, 'XTick', Np_vector(1:2:end))
