% plot simulation results
% y: NMSE
% x: SNR

close all;
%load results_11032015.mat
% load results_QPSK_training_11032015.mat
% load results_QPSK_Nb_2048_07082016.mat
% load results_Chu_Nb_2048_07082016.mat
% load results_QPSK_07082016.mat
load results_Chu_Np_2048_08162016.mat
load results_VAMP_Chu_Np_2048_01012017.mat
% load results_Chu_10db_Nb_2048_bit_1_inf_07112016.mat;
Nd = 16;
Nt = 64;
Nr = 64;

% Gnorm = sqrt(Nt*Nr) *ones(length(SNRdB_vector), length(bit_vector), length(Nb_vector), 10);
% NMSE_GAMP_dB = 10*log10(mean(MSE_GAMP, 4)./mean(Gnorm.^2, 4));

% NMSE_GAMP_dB = 10*log10(mean(MSE_GAMP, 4)/mean(Gnorm.^2));

% MSE_GAMP_dB = 10*log10(MSE_GAMP(:,:,:,4)./Gnorm(:,:,:,end));
% GI_MSE_GAMP_dB = 10*log10(mean(GI_MSE_GAMP, 4));

% GI_MSE_LS_dB = 10*log10(mean(GI_MSE_LS, 4));

% GI_MSE_LMMSE_dB = 10*log10(mean(GI_MSE_LMMSE, 4));

Training_length_index = 1;
% SNR_index = 1;

%% BG

% MSE versus SNR
figure,
plot1 = plot(SNRdB_vector, NMSE_BG_dB(:,:, Training_length_index), 'linewidth', 1);
set(plot1(1),'DisplayName','1-bit','Marker','o','LineStyle','-');
set(plot1(2),'DisplayName','2-bit','Marker','square','LineStyle','-');
set(plot1(3),'DisplayName','3-bit','Marker','diamond','LineStyle','-');
set(plot1(4),'DisplayName','4-bit','Marker','v','LineStyle','-');
set(plot1(5),'DisplayName','\infty-bit','Marker','^','LineStyle','--');
xlabel('SNR (dB)', 'fontsize',14)
ylabel('Normalized MSE (dB)', 'fontsize',14)
title('EM-BG-AMP algorithm')
grid on;
% legend('1-bit', '2-bit', '3-bit', '4-bit', '\infty-bit')

% % MSE versus Nb
% figure,
% plot1 = semilogx(Nb_vector, ( reshape( NMSE_GAMP_dB(SNR_index,:,:), length(bit_vector), []) )', 'linewidth', 1.5);
% set(plot1(1),'DisplayName','1-bit','Marker','o','LineStyle','-');
% set(plot1(2),'DisplayName','2-bit','Marker','square','LineStyle','-');
% set(plot1(3),'DisplayName','3-bit','Marker','diamond','LineStyle','-');
% set(plot1(4),'DisplayName','4-bit','Marker','v','LineStyle','-');
% set(plot1(5),'DisplayName','5-bit','Marker','^','LineStyle','--');
% xlabel('Training length', 'fontsize',14)
% ylabel('Normalized MSE (dB)', 'fontsize',14)
% %title('GAMP algorithm')
% grid on;
% legend('1 bit', '2 bits', '3 bits', '4 bits', '\infty bits')
% xlim([0 5000])
% set(gca, 'XTick', 0:1000:5000)



%% GM
% MSE versus SNR
figure,
plot1 = plot(SNRdB_vector, NMSE_GM_dB(:,:, Training_length_index), 'linewidth', 1);
set(plot1(1),'DisplayName','1-bit','Marker','o','LineStyle','-');
set(plot1(2),'DisplayName','2-bit','Marker','square','LineStyle','-');
set(plot1(3),'DisplayName','3-bit','Marker','diamond','LineStyle','-');
set(plot1(4),'DisplayName','4-bit','Marker','v','LineStyle','-');
set(plot1(5),'DisplayName','\infty-bit','Marker','^','LineStyle','--');
xlabel('SNR (dB)', 'fontsize',14)
ylabel('Normalized MSE (dB)', 'fontsize',14)
title('EM-GM-AMP algorithm')
grid on;

% % MSE versus Nb
% figure,
% plot1 = semilogx(Nb_vector, ( reshape( NMSE_GM_dB(SNR_index,:,:), length(bit_vector), []) )', 'linewidth', 1.5);
% set(plot1(1),'DisplayName','1-bit','Marker','o','LineStyle','-');
% set(plot1(2),'DisplayName','2-bit','Marker','square','LineStyle','-');
% set(plot1(3),'DisplayName','3-bit','Marker','diamond','LineStyle','-');
% set(plot1(4),'DisplayName','4-bit','Marker','v','LineStyle','-');
% set(plot1(5),'DisplayName','5-bit','Marker','^','LineStyle','--');
% xlabel('Training length', 'fontsize',14)
% ylabel('Normalized MSE (dB)', 'fontsize',14)
% %title('GAMP algorithm')
% grid on;
% legend('1 bit', '2 bits', '3 bits', '4 bits', '\infty bits')
% xlim([0 5000])
% set(gca, 'XTick', 0:1000:5000)


%% GAMP and VAMP

% MSE versus SNR
figure,
plot1 = plot(SNRdB_vector, NMSE_BG_dB(:,:, Training_length_index), 'linewidth', 1,'color','blue');
set(plot1(1),'DisplayName','1-bit','Marker','o','LineStyle','--');
set(plot1(2),'DisplayName','2-bit','Marker','square','LineStyle','--');
set(plot1(3),'DisplayName','3-bit','Marker','diamond','LineStyle','--');
set(plot1(4),'DisplayName','4-bit','Marker','v','LineStyle','--');
set(plot1(5),'DisplayName','\infty-bit','Marker','^','LineStyle','--');
xlabel('SNR (dB)', 'fontsize',14)
ylabel('Normalized MSE (dB)', 'fontsize',14)
grid on;
hold on;
plot1 = plot(SNRdB_vector, NMSE_GM_dB(:,:, Training_length_index), 'linewidth', 1,'color','black');
set(plot1(1),'DisplayName','1-bit','Marker','o','LineStyle','-');
set(plot1(2),'DisplayName','2-bit','Marker','square','LineStyle','-');
set(plot1(3),'DisplayName','3-bit','Marker','diamond','LineStyle','-');
set(plot1(4),'DisplayName','4-bit','Marker','v','LineStyle','-');
set(plot1(5),'DisplayName','\infty-bit','Marker','^','LineStyle','-');



plot1 = plot(SNRdB_vector, NMSE_VAMP_dB(:,:, Training_length_index), 'linewidth', 1,'color','red');
set(plot1(1),'DisplayName','1-bit','Marker','o','LineStyle','-');
set(plot1(2),'DisplayName','2-bit','Marker','square','LineStyle','-');
set(plot1(3),'DisplayName','3-bit','Marker','diamond','LineStyle','-');
set(plot1(4),'DisplayName','4-bit','Marker','v','LineStyle','-');
set(plot1(5),'DisplayName','\infty-bit','Marker','^','LineStyle','-');
legend('EM-BG-GAMP, 1-bit', 'EM-BG-GAMP, 2-bit','EM-BG-GAMP, 3-bit', 'EM-BG-GAMP, 4-bit','EM-BG-GAMP, \infty-bit', ...
    'EM-GM-GAMP, 1-bit','EM-GM-GAMP, 2-bit', 'EM-GM-GAMP, 3-bit','EM-GM-GAMP, 4-bit', 'EM-GM-GAMP, \infty-bit', ...
    'EM-GM-VAMP, 1-bit', 'EM-GM-VAMP, 2-bit','EM-GM-VAMP, 3-bit', 'EM-GM-VAMP, 4-bit','EM-GM-VAMP, \infty-bit',...
    'Location','best')


%% LS
% MSE versus SNR
figure,
plot1 = plot(SNRdB_vector, NMSE_LS_dB(:,:, Training_length_index), 'linewidth', 1);
set(plot1(1),'DisplayName','1-bit','Marker','o','LineStyle','-');
set(plot1(2),'DisplayName','2-bit','Marker','square','LineStyle','-');
set(plot1(3),'DisplayName','3-bit','Marker','diamond','LineStyle','-');
set(plot1(4),'DisplayName','4-bit','Marker','v','LineStyle','-');
set(plot1(5),'DisplayName','\infty-bit','Marker','^','LineStyle','--');
xlabel('SNR (dB)', 'fontsize',14)
ylabel('Normalized MSE (dB)', 'fontsize',14)
title('LS algorithm')
grid on;

%% LMMSE
% MSE versus SNR
figure,
plot1 = plot(SNRdB_vector, NMSE_LMMSE_dB(:,:, Training_length_index), 'linewidth', 1);
set(plot1(1),'DisplayName','1-bit','Marker','o','LineStyle','-');
set(plot1(2),'DisplayName','2-bit','Marker','square','LineStyle','-');
set(plot1(3),'DisplayName','3-bit','Marker','diamond','LineStyle','-');
set(plot1(4),'DisplayName','4-bit','Marker','v','LineStyle','-');
set(plot1(5),'DisplayName','\infty-bit','Marker','^','LineStyle','--');
xlabel('SNR (dB)', 'fontsize',14)
ylabel('Normalized MSE (dB)', 'fontsize',14)
title('ALMMSE algorithm')
grid on;


%% BPDN
% MSE versus SNR
figure,
plot1 = plot(SNRdB_vector, NMSE_BPDN_dB(:,:, Training_length_index), 'linewidth', 1);
set(plot1(1),'DisplayName','1-bit','Marker','o','LineStyle','-');
set(plot1(2),'DisplayName','2-bit','Marker','square','LineStyle','-');
set(plot1(3),'DisplayName','3-bit','Marker','diamond','LineStyle','-');
set(plot1(4),'DisplayName','4-bit','Marker','v','LineStyle','-');
set(plot1(5),'DisplayName','\infty-bit','Marker','^','LineStyle','--');
xlabel('SNR (dB)', 'fontsize',14)
ylabel('Normalized MSE (dB)', 'fontsize',14)
title('BPDN algorithm')
grid on;

%% QIHT
% MSE versus SNR
figure,
plot1 = plot(SNRdB_vector, NMSE_QIHT_dB(:,:, Training_length_index), 'linewidth', 1);
set(plot1(1),'DisplayName','1-bit','Marker','o','LineStyle','-');
set(plot1(2),'DisplayName','2-bit','Marker','square','LineStyle','-');
set(plot1(3),'DisplayName','3-bit','Marker','diamond','LineStyle','-');
set(plot1(4),'DisplayName','4-bit','Marker','v','LineStyle','-');
set(plot1(5),'DisplayName','\infty-bit','Marker','^','LineStyle','--');
xlabel('SNR (dB)', 'fontsize',14)
ylabel('Normalized MSE (dB)', 'fontsize',14)
title('QIHT algorithm')
grid on;

%% QIST
% MSE versus SNR
figure,
plot1 = plot(SNRdB_vector, NMSE_QIST_dB(:,:, Training_length_index), 'linewidth', 1);
set(plot1(1),'DisplayName','1-bit','Marker','o','LineStyle','-');
set(plot1(2),'DisplayName','2-bit','Marker','square','LineStyle','-');
set(plot1(3),'DisplayName','3-bit','Marker','diamond','LineStyle','-');
set(plot1(4),'DisplayName','4-bit','Marker','v','LineStyle','-');
set(plot1(5),'DisplayName','\infty-bit','Marker','^','LineStyle','--');
xlabel('SNR (dB)', 'fontsize',14)
ylabel('Normalized MSE (dB)', 'fontsize',14)
title('QIST algorithm')
grid on;

%% plot all the results in a single figure
figure, 
plot_LS = plot(SNRdB_vector, NMSE_LS_dB(:,:, Training_length_index), 'linewidth', 1,'Marker','d','color','[0 0.447 0.741]');
hold on;
grid on;
plot(SNRdB_vector, NMSE_LMMSE_dB(:,:, Training_length_index), 'linewidth', 1,'Marker','x','color','[0.85 0.325 0.098]');
plot(SNRdB_vector, NMSE_BPDN_dB(:,:, Training_length_index), 'linewidth', 1,'Marker','v','color','[0.929 0.694 0.125]');
plot(SNRdB_vector, NMSE_QIHT_dB(:,:, Training_length_index), 'linewidth', 1,'Marker','<','color','[0.494 0.184 0.556]');
plot(SNRdB_vector, NMSE_QIST_dB(:,:, Training_length_index), 'linewidth', 1,'Marker','^','color','[0.466 0.674 0.188]');
plot(SNRdB_vector, NMSE_BG_dB(:,:, Training_length_index), 'linewidth', 1,'Marker','o','color','[0.301 0.745 0.933]');
plot(SNRdB_vector, NMSE_GM_dB(:,:, Training_length_index), 'linewidth', 1,'Marker','s','color','[0.635 0.078 0.184]');
xlabel('SNR (dB)', 'fontsize',14)
ylabel('Normalized MSE (dB)', 'fontsize',14)




%% plot 1-bit results in a single figure
figure, 
plot_LS = plot(SNRdB_vector, NMSE_LS_dB(:,[1], Training_length_index), 'linewidth', 1,'Marker','d','color','[0 0.447 0.741]');
hold on;
grid on;
plot(SNRdB_vector, NMSE_LMMSE_dB(:,[1], Training_length_index), 'linewidth', 1,'Marker','x','color','[0.85 0.325 0.098]');
plot(SNRdB_vector, NMSE_BPDN_dB(:,[1], Training_length_index), 'linewidth', 1,'Marker','v','color','[0.929 0.694 0.125]');
plot(SNRdB_vector, NMSE_QIHT_dB(:,[1], Training_length_index), 'linewidth', 1,'Marker','<','color','[0.494 0.184 0.556]');
plot(SNRdB_vector, NMSE_QIST_dB(:,[1], Training_length_index), 'linewidth', 1,'Marker','^','color','[0.466 0.674 0.188]');
plot(SNRdB_vector, NMSE_BG_dB(:,[1], Training_length_index), 'linewidth', 1,'Marker','o','color','[0.301 0.745 0.933]');
plot(SNRdB_vector, NMSE_GM_dB(:,[1], Training_length_index), 'linewidth', 1,'Marker','s','color','[0.635 0.078 0.184]');

xlabel('SNR (dB)', 'fontsize',14)
ylabel('Normalized MSE (dB)', 'fontsize',14)
legend('LS, 1-bit','ALMMSE, 1-bit', 'BPDN, 1-bit', 'QIHT, 1-bit','QIST, 1-bit','BG, 1-bit', ...
    'GM, 1-bit', 'Location','best')


%% plot 4-bit results in a single figure
figure, 
plot_LS = plot(SNRdB_vector, NMSE_LS_dB(:,[4], Training_length_index), 'linewidth', 1,'Marker','d','color','[0 0.447 0.741]');
hold on;
grid on;
plot(SNRdB_vector, NMSE_LMMSE_dB(:,[4], Training_length_index), 'linewidth', 1,'Marker','x','color','[0.85 0.325 0.098]');
plot(SNRdB_vector, NMSE_BPDN_dB(:,[4], Training_length_index), 'linewidth', 1,'Marker','v','color','[0.929 0.694 0.125]');
plot(SNRdB_vector, NMSE_QIHT_dB(:,[4], Training_length_index), 'linewidth', 1,'Marker','<','color','[0.494 0.184 0.556]');
plot(SNRdB_vector, NMSE_BG_dB(:,[4], Training_length_index), 'linewidth', 1,'Marker','o','color','[0.301 0.745 0.933]');
plot(SNRdB_vector, NMSE_GM_dB(:,[4], Training_length_index), 'linewidth', 1,'Marker','s','color','[0.635 0.078 0.184]');
plot(SNRdB_vector, NMSE_VAMP_dB(:,[4], Training_length_index), 'linewidth', 1,'Marker','^','color','[0.466 0.674 0.188]');
xlabel('SNR (dB)', 'fontsize',14)
ylabel('Normalized MSE (dB)', 'fontsize',14)
legend('LS','ALMMSE', 'SPGL1', 'QIHT','EM-BG-GAMP', ...
    'EM-GM-GAMP', 'EM-GM-VAMP', 'Location','best')


%% plot 4-bit and infinite-bit results in a single figure
figure, 
plot_LS = plot(SNRdB_vector, NMSE_LS_dB(:,[4, 5], Training_length_index), 'linewidth', 1,'Marker','d','color','[0 0.447 0.741]');
hold on;
grid on;
plot(SNRdB_vector, NMSE_LMMSE_dB(:,[4, 5], Training_length_index), 'linewidth', 1,'Marker','x','color','[0.85 0.325 0.098]');
plot(SNRdB_vector, NMSE_BPDN_dB(:,[4, 5], Training_length_index), 'linewidth', 1,'Marker','v','color','[0.929 0.694 0.125]');
plot(SNRdB_vector, NMSE_QIHT_dB(:,[4, 5], Training_length_index), 'linewidth', 1,'Marker','<','color','[0.494 0.184 0.556]');
plot(SNRdB_vector, NMSE_BG_dB(:,[4, 5], Training_length_index), 'linewidth', 1,'Marker','o','color','[0.301 0.745 0.933]');
plot(SNRdB_vector, NMSE_GM_dB(:,[4, 5], Training_length_index), 'linewidth', 1,'Marker','s','color','[0.635 0.078 0.184]');
plot(SNRdB_vector, NMSE_VAMP_dB(:,[4, 5], Training_length_index), 'linewidth', 1,'Marker','^','color','[0.466 0.674 0.188]');

xlabel('SNR (dB)', 'fontsize',14)
ylabel('Normalized MSE (dB)', 'fontsize',14)
legend('LS, 4-bit', 'LS, \infty-bit','ALMMSE, 4-bit', 'ALMMSE, \infty-bit','SPGL1, 4-bit', 'SPGL1, \infty-bit',...
    'QIHT, 4-bit','QIHT, \infty-bit', 'EM-BG-GAMP, 4-bit', 'EM-BG-GAMP, \infty-bit',...
    'EM-GM-GAMP, 4-bit', 'EM-GM-GAMP, \infty-bit', 'EM-GM-VAMP, 4-bit','EM-GM-VAMP, \infty-bit', 'Location','best')




load results_GAMP_VAMP_Chu_Np_2048_02252017.mat

NMSE_BG_VAMP_dB = Results.NMSE.NMSE_BG_VAMP_dB;
NMSE_GM_VAMP_dB = Results.NMSE.NMSE_GM_VAMP_dB;
NMSE_BG_GAMP_dB = Results.NMSE.NMSE_BG_GAMP_dB;
NMSE_GM_GAMP_dB = Results.NMSE.NMSE_GM_GAMP_dB;

SNRdB_vector = Setup.SNRdB_vector;

% MSE versus SNR
figure,
plot1 = plot(SNRdB_vector, NMSE_BG_GAMP_dB(:,:, Training_length_index), 'linewidth', 1,'color','m');
set(plot1(1),'DisplayName','1-bit','Marker','o','LineStyle','--');
set(plot1(2),'DisplayName','2-bit','Marker','square','LineStyle','--');
set(plot1(3),'DisplayName','3-bit','Marker','diamond','LineStyle','--');
set(plot1(4),'DisplayName','4-bit','Marker','v','LineStyle','--');
set(plot1(5),'DisplayName','\infty-bit','Marker','^','LineStyle','--');
xlabel('SNR (dB)', 'fontsize',14)
ylabel('Normalized MSE (dB)', 'fontsize',14)
grid on;
hold on;
plot1 = plot(SNRdB_vector, NMSE_GM_GAMP_dB(:,:, Training_length_index), 'linewidth', 1,'color','red');
set(plot1(1),'DisplayName','1-bit','Marker','o','LineStyle','-');
set(plot1(2),'DisplayName','2-bit','Marker','square','LineStyle','-');
set(plot1(3),'DisplayName','3-bit','Marker','diamond','LineStyle','-');
set(plot1(4),'DisplayName','4-bit','Marker','v','LineStyle','-');
set(plot1(5),'DisplayName','\infty-bit','Marker','^','LineStyle','-');

plot1 = plot(SNRdB_vector, NMSE_BG_VAMP_dB(:,:, Training_length_index), 'linewidth', 1,'color','blue');
set(plot1(1),'DisplayName','1-bit','Marker','o','LineStyle','--');
set(plot1(2),'DisplayName','2-bit','Marker','square','LineStyle','--');
set(plot1(3),'DisplayName','3-bit','Marker','diamond','LineStyle','--');
set(plot1(4),'DisplayName','4-bit','Marker','v','LineStyle','--');
set(plot1(5),'DisplayName','\infty-bit','Marker','^','LineStyle','--');

plot1 = plot(SNRdB_vector, NMSE_GM_VAMP_dB(:,:, Training_length_index), 'linewidth', 1,'color','black');
set(plot1(1),'DisplayName','1-bit','Marker','o','LineStyle','-');
set(plot1(2),'DisplayName','2-bit','Marker','square','LineStyle','-');
set(plot1(3),'DisplayName','3-bit','Marker','diamond','LineStyle','-');
set(plot1(4),'DisplayName','4-bit','Marker','v','LineStyle','-');
set(plot1(5),'DisplayName','\infty-bit','Marker','^','LineStyle','-');
legend({'EM-BG-GAMP, 1-bit', 'EM-BG-GAMP, 2-bit','EM-BG-GAMP, 3-bit', 'EM-BG-GAMP, 4-bit','EM-BG-GAMP, \infty-bit', ...
    'EM-GM-GAMP, 1-bit','EM-GM-GAMP, 2-bit', 'EM-GM-GAMP, 3-bit','EM-GM-GAMP, 4-bit', 'EM-GM-GAMP, \infty-bit', ...
    'EM-BG-VAMP, 1-bit', 'EM-BG-VAMP, 2-bit','EM-BG-VAMP, 3-bit', 'EM-BG-VAMP, 4-bit','EM-BG-VAMP, \infty-bit', ...
    'EM-GM-VAMP, 1-bit', 'EM-GM-VAMP, 2-bit','EM-GM-VAMP, 3-bit', 'EM-GM-VAMP, 4-bit','EM-GM-VAMP, \infty-bit'},...
    'Location','best')
