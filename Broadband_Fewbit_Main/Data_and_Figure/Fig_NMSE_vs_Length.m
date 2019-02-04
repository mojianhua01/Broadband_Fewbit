% plot simulation results
% y: NMSE
% x: Training Length

close all;
%load results_11032015.mat
% load results_QPSK_training_11032015.mat
% load results_QPSK_Nb_2048_07082016.mat
% load results_Chu_Nb_2048_07082016.mat
% load results_QPSK_07082016.mat
% load results_Chu_0dB_Nb_512_4096_07112016.mat
load results_Chu_0dB_Np_1024_5120_08212016.mat
% load results_Chu_10dB_Np_1024_5120_09232016.mat
% load results_VAMP_Chu_0dB_Np_1024_5120_01012017.mat
load results_VAMP_Chu_0dB_Np_1024_5120_01012017.mat
Nd = 16;
Nt = 64;
Nr = 64;

% Training_length_index = 1;
SNR_index = 1;

% bit_vector = [1 2 3 4];

%% BG

% % MSE versus Np
figure,
plot1 = semilogx(Np_vector, ( reshape( NMSE_BG_dB(SNR_index,:,:), length(bit_vector), []) )', 'linewidth', 1.5);
set(plot1(1),'DisplayName','1-bit','Marker','o','LineStyle','-');
set(plot1(2),'DisplayName','2-bit','Marker','square','LineStyle','-');
set(plot1(3),'DisplayName','3-bit','Marker','diamond','LineStyle','-');
set(plot1(4),'DisplayName','4-bit','Marker','v','LineStyle','-');
set(plot1(5),'DisplayName','\infty-bit','Marker','^','LineStyle','--');
xlabel('Training Length', 'fontsize',14)
ylabel('Normalized MSE (dB)', 'fontsize',14)
title('EM-BG-AMP algorithm')
grid on;
xlim([0 5000])
set(gca, 'XTick', Np_vector)


%% GM

% % MSE versus Np
figure,
plot1 = semilogx(Np_vector, ( reshape( NMSE_GM_dB(SNR_index,:,:), length(bit_vector), []) )', 'linewidth', 1.5);
set(plot1(1),'DisplayName','1-bit','Marker','o','LineStyle','-');
set(plot1(2),'DisplayName','2-bit','Marker','square','LineStyle','-');
set(plot1(3),'DisplayName','3-bit','Marker','diamond','LineStyle','-');
set(plot1(4),'DisplayName','4-bit','Marker','v','LineStyle','-');
set(plot1(5),'DisplayName','\infty-bit','Marker','^','LineStyle','--');
xlabel('Training Length', 'fontsize',14)
ylabel('Normalized MSE (dB)', 'fontsize',14)
title('EM-GM-AMP algorithm')
grid on;
xlim([0 5000])
set(gca, 'XTick', Np_vector)


% %% BPDN
% 
% % % MSE versus Np
% figure,
% plot1 = semilogx(Np_vector, ( reshape( NMSE_BPDN_dB(SNR_index,:,:), length(bit_vector), []) )', 'linewidth', 1.5);
% set(plot1(1),'DisplayName','1-bit','Marker','o','LineStyle','-');
% set(plot1(2),'DisplayName','2-bit','Marker','square','LineStyle','-');
% set(plot1(3),'DisplayName','3-bit','Marker','diamond','LineStyle','-');
% set(plot1(4),'DisplayName','4-bit','Marker','v','LineStyle','-');
% set(plot1(5),'DisplayName','\infty-bit','Marker','^','LineStyle','--');
% xlabel('Training Length', 'fontsize',14)
% ylabel('Normalized MSE (dB)', 'fontsize',14)
% title('BPDN algorithm')
% grid on;
% xlim([0 5000])
% set(gca, 'XTick', Np_vector)
% 
% %% LS
% 
% % % MSE versus Np
% figure,
% plot1 = semilogx(Np_vector, ( reshape( NMSE_LS_dB(SNR_index,:,:), length(bit_vector), []) )', 'linewidth', 1.5);
% set(plot1(1),'DisplayName','1-bit','Marker','o','LineStyle','-');
% set(plot1(2),'DisplayName','2-bit','Marker','square','LineStyle','-');
% set(plot1(3),'DisplayName','3-bit','Marker','diamond','LineStyle','-');
% set(plot1(4),'DisplayName','4-bit','Marker','v','LineStyle','-');
% set(plot1(5),'DisplayName','\infty-bit','Marker','^','LineStyle','--');
% xlabel('Training Length', 'fontsize',14)
% ylabel('Normalized MSE (dB)', 'fontsize',14)
% title('LS algorithm')
% grid on;
% xlim([0 5000])
% set(gca, 'XTick', Np_vector)
% 
% 
% %% LMMSE
% 
% % % MSE versus Np
% figure,
% plot1 = semilogx(Np_vector, ( reshape( NMSE_LMMSE_dB(SNR_index,:,:), length(bit_vector), []) )', 'linewidth', 1.5);
% set(plot1(1),'DisplayName','1-bit','Marker','o','LineStyle','-');
% set(plot1(2),'DisplayName','2-bit','Marker','square','LineStyle','-');
% set(plot1(3),'DisplayName','3-bit','Marker','diamond','LineStyle','-');
% set(plot1(4),'DisplayName','4-bit','Marker','v','LineStyle','-');
% set(plot1(5),'DisplayName','\infty-bit','Marker','^','LineStyle','--');
% xlabel('Training Length', 'fontsize',14)
% ylabel('Normalized MSE (dB)', 'fontsize',14)
% title('LMMSE algorithm')
% grid on;
% xlim([0 5000])
% set(gca, 'XTick', Np_vector)


%% plot 1-bit and infinite-bit results in a single figure
figure,
plot1 = semilogx(Np_vector, ( reshape( NMSE_BG_dB(SNR_index,:,:), length(bit_vector), []) )', 'linewidth', 1, 'color', 'blue');
set(plot1(1),'DisplayName','1-bit','Marker','o','LineStyle','--');
set(plot1(2),'DisplayName','2-bit','Marker','square','LineStyle','--');
set(plot1(3),'DisplayName','3-bit','Marker','diamond','LineStyle','--');
set(plot1(4),'DisplayName','4-bit','Marker','v','LineStyle','--');
set(plot1(5),'DisplayName','\infty-bit','Marker','^','LineStyle','--');
hold on;
grid on;
plot1 = semilogx(Np_vector, ( reshape( NMSE_GM_dB(SNR_index,:,:), length(bit_vector), []) )', 'linewidth', 1, 'color', 'black');
set(plot1(1),'DisplayName','1-bit','Marker','o','LineStyle','-');
set(plot1(2),'DisplayName','2-bit','Marker','square','LineStyle','-');
set(plot1(3),'DisplayName','3-bit','Marker','diamond','LineStyle','-');
set(plot1(4),'DisplayName','4-bit','Marker','v','LineStyle','-');
set(plot1(5),'DisplayName','\infty-bit','Marker','^','LineStyle','-');

plot1 = semilogx(Np_vector, ( reshape( NMSE_VAMP_dB(SNR_index,:,:), length(bit_vector), []) )', 'linewidth', 1, 'color', 'red');
set(plot1(1),'DisplayName','1-bit','Marker','o','LineStyle','-');
set(plot1(2),'DisplayName','2-bit','Marker','square','LineStyle','-');
set(plot1(3),'DisplayName','3-bit','Marker','diamond','LineStyle','-');
set(plot1(4),'DisplayName','4-bit','Marker','v','LineStyle','-');
set(plot1(5),'DisplayName','\infty-bit','Marker','^','LineStyle','-');

xlabel('Training Length', 'fontsize',14)
ylabel('Normalized MSE (dB)', 'fontsize',14)
legend('EM-BG-GAMP, 1-bit', 'EM-BG-GAMP, 2-bit','EM-BG-GAMP, 3-bit', 'EM-BG-GAMP, 4-bit','EM-BG-GAMP, \infty-bit', ...
    'EM-GM-GAMP, 1-bit','EM-GM-GAMP, 2-bit', 'EM-GM-GAMP, 3-bit','EM-GM-GAMP, 4-bit', 'EM-GM-GAMP, \infty-bit', ...
    'EM-GM-VAMP, 1-bit', 'EM-GM-VAMP, 2-bit','EM-GM-VAMP, 3-bit', 'EM-GM-VAMP, 4-bit','EM-GM-VAMP, \infty-bit',...
    'Location','best')
xlim([0 9000])
set(gca, 'XTick', Np_vector(1:2:end))



load results_Chu_0dB_VAMP_GAMP_Np_1024_5120_02252017.mat
SNR_index = 1;
Rate_BG_VAMP_mean = Results.Rate.Rate_BG_VAMP_mean;
Rate_GM_VAMP_mean = Results.Rate.Rate_GM_VAMP_mean;
Rate_BG_GAMP_mean = Results.Rate.Rate_BG_GAMP_mean;
Rate_GM_GAMP_mean = Results.Rate.Rate_GM_GAMP_mean;

NMSE_BG_VAMP_dB = Results.NMSE.NMSE_BG_VAMP_dB;
NMSE_GM_VAMP_dB = Results.NMSE.NMSE_GM_VAMP_dB;
NMSE_BG_GAMP_dB = Results.NMSE.NMSE_BG_GAMP_dB;
NMSE_GM_GAMP_dB = Results.NMSE.NMSE_GM_GAMP_dB;

Np_vector = Setup.Np_vector;
bit_vector = Setup.bit_vector;
%% plot 1-4 bit and infinite-bit results in a single figure
figure,
plot1 = semilogx(Np_vector, ( reshape( NMSE_BG_GAMP_dB(SNR_index,:,:), length(bit_vector), []) )', 'linewidth', 1, 'color', 'magenta');
set(plot1(1),'DisplayName','1-bit','Marker','o','LineStyle','--');
set(plot1(2),'DisplayName','2-bit','Marker','square','LineStyle','--');
set(plot1(3),'DisplayName','3-bit','Marker','diamond','LineStyle','--');
set(plot1(4),'DisplayName','4-bit','Marker','v','LineStyle','--');
set(plot1(5),'DisplayName','\infty-bit','Marker','^','LineStyle','--');
hold on;
grid on;
plot1 = semilogx(Np_vector, ( reshape( NMSE_GM_GAMP_dB(SNR_index,:,:), length(bit_vector), []) )', 'linewidth', 1, 'color', 'red');
set(plot1(1),'DisplayName','1-bit','Marker','o','LineStyle','-');
set(plot1(2),'DisplayName','2-bit','Marker','square','LineStyle','-');
set(plot1(3),'DisplayName','3-bit','Marker','diamond','LineStyle','-');
set(plot1(4),'DisplayName','4-bit','Marker','v','LineStyle','-');
set(plot1(5),'DisplayName','\infty-bit','Marker','^','LineStyle','-');

plot1 = semilogx(Np_vector, ( reshape( NMSE_BG_VAMP_dB(SNR_index,:,:), length(bit_vector), []) )', 'linewidth', 1, 'color', 'blue');
set(plot1(1),'DisplayName','1-bit','Marker','o','LineStyle','--');
set(plot1(2),'DisplayName','2-bit','Marker','square','LineStyle','--');
set(plot1(3),'DisplayName','3-bit','Marker','diamond','LineStyle','--');
set(plot1(4),'DisplayName','4-bit','Marker','v','LineStyle','--');
set(plot1(5),'DisplayName','\infty-bit','Marker','^','LineStyle','--');

plot1 = semilogx(Np_vector, ( reshape( NMSE_GM_VAMP_dB(SNR_index,:,:), length(bit_vector), []) )', 'linewidth', 1, 'color', 'black');
set(plot1(1),'DisplayName','1-bit','Marker','o','LineStyle','-');
set(plot1(2),'DisplayName','2-bit','Marker','square','LineStyle','-');
set(plot1(3),'DisplayName','3-bit','Marker','diamond','LineStyle','-');
set(plot1(4),'DisplayName','4-bit','Marker','v','LineStyle','-');
set(plot1(5),'DisplayName','\infty-bit','Marker','^','LineStyle','-');


xlabel('Training Length', 'fontsize',14)
ylabel('Normalized MSE (dB)', 'fontsize',14)
legend({'EM-BG-GAMP, 1-bit', 'EM-BG-GAMP, 2-bit','EM-BG-GAMP, 3-bit', 'EM-BG-GAMP, 4-bit','EM-BG-GAMP, \infty-bit', ...
    'EM-GM-GAMP, 1-bit','EM-GM-GAMP, 2-bit', 'EM-GM-GAMP, 3-bit','EM-GM-GAMP, 4-bit', 'EM-GM-GAMP, \infty-bit', ...
    'EM-BG-VAMP, 1-bit', 'EM-BG-VAMP, 2-bit','EM-BG-VAMP, 3-bit', 'EM-BG-VAMP, 4-bit','EM-BG-VAMP, \infty-bit', ...
    'EM-GM-VAMP, 1-bit', 'EM-GM-VAMP, 2-bit','EM-GM-VAMP, 3-bit', 'EM-GM-VAMP, 4-bit','EM-GM-VAMP, \infty-bit'},...
    'Location','best')
xlim([0 9000])
set(gca, 'XTick', Np_vector(1:2:end))