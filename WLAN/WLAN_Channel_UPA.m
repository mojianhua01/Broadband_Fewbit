% function [H_scaled, H_v] = func_mmWave_Channel_UPA(BS_ant_W, BS_ant_H, MS_ant_W,MS_ant_H, ...
%     Angle_spread_W, Angle_spread_H, Num_cluster, Num_ray, Narrowband_flag, Nd)

%This code generates the 60GHz WLAN channel

%% setup for test
clear;
clc;
close all;

BS_ant_W = 4;
BS_ant_H = 4;
MS_ant_W = 4;
MS_ant_H = 4;

% BS_ant_W = 1;
% BS_ant_H = 1;
% MS_ant_W = 1;
% MS_ant_H = 1;

% Simulation parameters
BS_ant = BS_ant_W * BS_ant_H;
MS_ant = MS_ant_W * MS_ant_H;
BS_ant_W_index=0:1:BS_ant_W-1;
BS_ant_H_index=0:1:BS_ant_H-1;
MS_ant_W_index=0:1:MS_ant_W-1;
MS_ant_H_index=0:1:MS_ant_H-1;
%     H=zeros(MS_ant,BS_ant);
KD=pi;  % Assuming: K=2pi/lambda, D=lambda/2
for channel_index=1:1:100
    
    % Channel generation
    % load configuration structure <- cr_ch_cfg.m
    cfg = cr_ch_cfg;
    
    % generate space-time channel impulse response realization
    ch = gen_cr_ch(cfg.cr,cfg.bf.ps,cfg.bf.pol);
    
    % The following code are used to generate a SISO channel
    % apply beamforming algorithm
    [imp_res,toa] = beamforming(cfg.bf,ch);
    
    % continuous time to descrete time conversion
    imp_res = ct2dt(imp_res,toa,cfg.sample_rate);
    
    % normalization according to Pnorm parameter
    if (cfg.Pnorm)
        imp_res = imp_res./norm(imp_res);
    end
    
%     figure,
%     bar(abs(imp_res));
%     title('WLAN MIMO Channel With Beamforming')
    
    Num_ray = length(ch.am);
    H_ct = zeros(MS_ant, BS_ant, length(ch.am));
    
    for ray_ix = 1:1:length(ch.am)
        alpha       = ch.am(ray_ix, 1); % gain
        AoD_phi     = ch.tx_az(ray_ix, 1);
        AoD_theta   = ch.tx_el(ray_ix, 1);
        AoA_phi     = ch.rx_az(ray_ix, 1);
        AoA_theta   = ch.rx_el(ray_ix, 1);
        
        BS_steering = sqrt(1/BS_ant)* transpose( exp(1j*KD*BS_ant_W_index*sin(AoD_phi)*sin(AoD_theta)) ) * exp(1j*KD*BS_ant_H_index*cos(AoD_theta));
        
        MS_steering = sqrt(1/MS_ant)* transpose( exp(1j*KD*MS_ant_W_index*sin(AoA_phi)*sin(AoA_theta)) ) * exp(1j*KD*MS_ant_H_index*cos(AoA_theta));
        
        H_ct(:,:, ray_ix) = alpha * MS_steering(:) * transpose( BS_steering(:) );
    end;
    
    t_s = 1./cfg.sample_rate;
    
    t_0 = ch.toa - min(ch.toa);
    
    Nd = round(max(t_0)./t_s) + 1;
    
    H_dt = zeros(MS_ant, BS_ant, Nd);
    
    for ray_ix = 1:length(t_0)
        
        time_bin = round(t_0(ray_ix)./t_s) + 1;
        
        H_dt(:, :, time_bin) = H_dt(:, :, time_bin) + H_ct(:,:, ray_ix);
    end
    
    Hv = zeros(MS_ant, BS_ant, Nd);
    
    for d=1:Nd
        Hv(:,:,d) = 1/sqrt(BS_ant * MS_ant) * dftmtx(MS_ant)'* H_dt(:,:,d) *dftmtx(BS_ant);
    end
    
    Hv = Hv./norm( reshape(Hv,[],1)) * sqrt(BS_ant * MS_ant);
    
    Hv_all{channel_index} = Hv;
    
end;
save H_WLAN_16_16.mat Hv_all;
if BS_ant ~= 1 || MS_ant ~= 1
    figure,
    bar3((abs(reshape(permute(Hv,[1,3,2]), MS_ant*Nd, BS_ant))))
    set(gca, 'FontSize',12)
    grid on;
    xlabel('Tx angle')
    ylabel('Rx angle & delay')
    zlabel('$|\mathbf{H}_{\mathrm{v}}|$', 'Interpreter', 'latex');
    title('WLAN MIMO Channel Without Beamforming')
end;

if BS_ant == 1 && MS_ant == 1
    figure,
    bar(abs(reshape(permute(Hv,[1,3,2]), 1, [])));
    title('WLAN SISO Channel')
end;


% run_flag = 1
%
% if run_flag == 0
%     return;
% end
%
% % Model parameters
%
% % Channel model
% Nt_W = BS_ant_W;
% Nt_H = BS_ant_H;
% Nt = Nt_W * Nt_H;
% Nr_W = MS_ant_W;
% Nr_H = MS_ant_H;
% Nr = Nr_W * Nr_H;
%
% Hybrid_flag = 0;
%
% % ADC setup
% bit_vector = 8; %[1, 2, 3, 4, 5, 6, 7, 8, +inf];
%
% % Pilot length and SNR
% N = 2^11;
% SNR_vector = 20;
%
% % EG-BG-GAMP setup
% BGmean = 0;        	% Bernoulli-Gaussian active mean
% BGvar  = 1;         	% Bernoulli-Gaussian active variance
%
% % Walsh_Hadamard_matrix = hadamard(Num_cluster);
% % S_train = Walsh_Hadamard_matrix(randsample(1:Num_cluster,Nt), 1:Num_cluster);
%
% % S_train = (randn(Nt, N) + 1j * randn(Nt, N));
% S_train = sign(randn(Nt,N))/sqrt(2)+1i*sign(randn(Nt,N))/sqrt(2);
%
% Ns = 4;
% if Hybrid_flag == 1
%     Analog_precoding_pattern = zeros(Nt, N);
%
%     for i=1:1:N
%         temp = randperm(Nt, Ns)';
%         Analog_precoding_pattern(temp, i) = 1;
%     end;
%
%     S_train_hybrid = Analog_precoding_pattern.*S_train;
%
%     for i=1:1:N
%         temp = randperm(Nr, Ns)';
%         Analog_combining_pattern(temp, i) = 1;
%     end;
% end;
%
% for index_channel = 1:1:1
%     %     Gtrue = func_mmWave_Channel_UPA(Nt_W, Nt_H, Nr_W, Nr_H, ...
%     %         Angle_spread_W, Angle_spread_H, Num_cluster, Num_ray, Narrowband_flag, Nd);
%     Gtrue = Hv_all{index_channel};
%     gtrue = reshape(Gtrue, [],1);
%
%     for index_bit = 1:1:length(bit_vector)
%
%         bit = bit_vector(index_bit);
%
%         for index_SNR = 1:1:length(SNR_vector)
%
%             SNR = SNR_vector(index_SNR);
%
%             if Hybrid_flag == 0
%                 [ Ghat ] = func_Broadband_Few_Bit( Gtrue, SNR, S_train, bit);
%             else
%                 [ Ghat ] = func_Broadband_Few_Bit_Hybrid (Gtrue, SNR, S_train_hybrid, bit, Analog_combining_pattern);
%             end;
%
%             ghat = reshape(Ghat, [], 1);
%
%             GAMP_MSE(index_SNR, index_bit, index_channel) = norm(gtrue-ghat, 2)^2;
%
%             Gnorm(index_SNR, index_bit, index_channel) = norm(gtrue)^2;
%         end;
%     end;
% end;
%
%
% figure,
% plot(SNR_vector, 10*log10(mean(GAMP_MSE, 3)./mean(Gnorm, 3)), 'r');
% xlabel('SNR', 'fontsize',14)
% ylabel('Normalized MSE', 'fontsize',14)
%
% 10*log10(mean(GAMP_MSE, 3)./mean(Gnorm, 3))
%
% figure,
% Nd = length(gtrue)/Nt/Nr;
% subplot(121)
% gtrue2 = reshape(permute(reshape(gtrue,Nr,Nt,Nd),[1,3,2]),Nd*Nr,Nt);
% bar3(abs(gtrue2));
% title('True value of the channel', 'Fontsize', 14)
% xlabel('Tx angle', 'Fontsize', 12)
% ylabel('Rx angle & delay', 'Fontsize', 12)
%
% subplot(122)
% gtrue2 = reshape(permute(reshape(ghat,Nr,Nt,Nd),[1,3,2]),Nd*Nr,Nt);
% bar3(abs(gtrue2));
% title('Estimate of the channel', 'Fontsize', 14)
% xlabel('Tx angle', 'Fontsize', 12)
% ylabel('Rx angle & delay', 'Fontsize', 12)
