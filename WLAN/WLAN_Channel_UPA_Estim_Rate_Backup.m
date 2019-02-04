clc;
clear all;
close all;
% Model parameters

% Channel model
Nt_W = 4;
Nt_H = 4;
Nt = Nt_W * Nt_H;
Nr_W = 4;
Nr_H = 4;
Nr = Nr_W * Nr_H;

Hybrid_flag = 1;

% ADC setup
bit_vector = [1, 2, 3, 4, 5, 6, 7, 8, +inf];

beta_vector = [0.3634 0.1175 0.03454 0.009497 0.002499 0.0006642 0.000166 0.00004151];

% Pilot length and SNR
N = 2^12;
SNR_vector = -10:2:10;

% EG-BG-GAMP setup
BGmean = 0;        	% Bernoulli-Gaussian active mean
BGvar  = 1;         	% Bernoulli-Gaussian active variance

% Walsh_Hadamard_matrix = hadamard(Num_cluster);
% S_train = Walsh_Hadamard_matrix(randsample(1:Num_cluster,Nt), 1:Num_cluster);

% S_train = (randn(Nt, N) + 1j * randn(Nt, N));
S_train = sign(randn(Nt,N))/sqrt(2)+1i*sign(randn(Nt,N))/sqrt(2);

Ns = 2;
if Hybrid_flag == 1
    Analog_precoding_pattern = zeros(Nt, N);
    
    for i=1:1:N
        temp = randperm(Nt, Ns)';
        Analog_precoding_pattern(temp, i) = 1;
    end;
    
    S_train_hybrid = Analog_precoding_pattern.*S_train;
    
    for i=1:1:N
        temp = randperm(Nr, Ns)';
        Analog_combining_pattern(temp, i) = 1;
    end;
end;

load H_WLAN_16_16.mat

for index_channel = 1:1:20
    %     Gtrue = func_mmWave_Channel_UPA(Nt_W, Nt_H, Nr_W, Nr_H, ...
    %         Angle_spread_W, Angle_spread_H, Num_cluster, Num_ray, Narrowband_flag, Nd);
    Gtrue = Hv_all{index_channel};
    gtrue = reshape(Gtrue, [],1);
    
    for index_bit = 1:1:length(bit_vector)
        
        bit = bit_vector(index_bit);
        
        for index_SNR = 1:1:length(SNR_vector)
            
            SNR = SNR_vector(index_SNR);
            
            if Hybrid_flag == 0
                [ Ghat_GAMP ]   = func_Broadband_Few_Bit( Gtrue, SNR, S_train, bit);
%                 [ Ghat_LS ]     = func_Broadband_Few_Bit_LS( Gtrue, SNR, S_train, bit);
            else
                [ Ghat_GAMP ]   = func_Broadband_Few_Bit_Hybrid (Gtrue, SNR, S_train_hybrid, bit, Analog_combining_pattern);
            end;
            ghat_GAMP = reshape(Ghat_GAMP, [], 1);
            MSE_GAMP(index_SNR, index_bit, index_channel) = norm(gtrue-ghat_GAMP, 2)^2;
            
%             ghat_LS = reshape(Ghat_LS, [], 1);
%             MSE_LS(index_SNR, index_bit, index_channel) = norm(gtrue-ghat_LS, 2)^2;
            
            Gnorm(index_SNR, index_bit, index_channel) = norm(gtrue)^2;
                
            Ghat_GAMP_Angular = sum( abs(Ghat_GAMP).^2, 3);
            
            Ghat_GAMP_Angular_reduced = Ghat_GAMP_Angular;
            
            for index_Ns = 1:1:Ns
                [AoA(index_Ns), AoD(index_Ns)] = ...
                    find(Ghat_GAMP_Angular == max(max(Ghat_GAMP_Angular_reduced)));
                Ghat_GAMP_Angular_reduced = Ghat_GAMP_Angular( setdiff([1:1:Nr], AoA), ...
                    setdiff([1:1:Nt], AoD));
            end;
            
            Ghat_GAMP_Effective = Ghat_GAMP(AoA,AoD, :) ;
            Gtrue_Effective = Gtrue(AoA, AoD, :) ;

            % ???????????????
            for index_Ns_r  = 1:1:Ns
                for index_Ns_t = 1:1:Ns
                    Ghat_GAMP_Effective_FD(index_Ns_r, index_Ns_t,:) = ...
                        fft( reshape(Ghat_GAMP_Effective(index_Ns_r, index_Ns_t, :),1,[]), 256);
                    Gtrue_Effective_FD(index_Ns_r, index_Ns_t,:) = ...
                        fft( reshape(Gtrue_Effective(index_Ns_r, index_Ns_t, :),1,[]), 256);
                end;
            end;
            
            Gerror_GAMP_Effective_FD = Ghat_GAMP_Effective_FD - Gtrue_Effective_FD;
            
            if bit == +inf
                beta = 0;
            else
                beta = beta_vector(bit);
            end;
            
            for sub_index = 1:1:256
                Rate_sub(sub_index) = log2( det ( eye(Ns) + ...
                    (1- beta ) * Ghat_GAMP_Effective_FD(:,:,sub_index) * ...
                    Ghat_GAMP_Effective_FD(:,:,sub_index)' * ...
                    ( (1-beta)* Gerror_GAMP_Effective_FD(:,:,sub_index) * Gerror_GAMP_Effective_FD(:,:,sub_index)' ...
                    + Ns/(10^(SNR/10)) + ...
                    beta * Gtrue_Effective_FD(:,:,sub_index)*Gtrue_Effective_FD(:,:,sub_index)' )^(-1)));
            end;
            Rate(index_SNR, index_bit, index_channel) = real(sum(Rate_sub)); 
            
        end;
    end;
end;

%%%% GAMP
Nd = size(Gtrue, 3);

figure,
plot(SNR_vector, 10*log10(mean(MSE_GAMP, 3)./mean(Gnorm, 3)), 'r');
xlabel('SNR', 'fontsize',14)
ylabel('Normalized MSE', 'fontsize',14)
title('GAMP algorithm')

10*log10(mean(MSE_GAMP, 3)./mean(Gnorm, 3))

figure,
subplot(121)
gtrue2 = reshape(permute(reshape(gtrue,Nr,Nt,Nd),[1,3,2]),Nd*Nr,Nt);
bar3(abs(gtrue2));
title('True value of the channel', 'Fontsize', 14)
xlabel('Tx angle', 'Fontsize', 12)
ylabel('Rx angle & delay', 'Fontsize', 12)

subplot(122)
gtrue2 = reshape(permute(reshape(ghat_GAMP,Nr,Nt,Nd),[1,3,2]),Nd*Nr,Nt);
bar3(abs(gtrue2));
title('Estimate of the channel', 'Fontsize', 14)
xlabel('Tx angle', 'Fontsize', 12)
ylabel('Rx angle & delay', 'Fontsize', 12)

%% Delete the bad value
ii = 1;
for index_channel = 1:1:20
    if isequal (sort( reshape(MSE_GAMP(:,:, index_channel), length(SNR_vector), []), 'descend' ) ...
            , reshape(MSE_GAMP(:,:, index_channel), length(SNR_vector), []))
        MSE_GAMP_M(:,:, ii) = MSE_GAMP(:,:, index_channel);
        ii = ii + 1;
    end;
end;

Outlier_prob_GAMP = (ii-1)/20

figure,
plot(SNR_vector, 10*log10(mean(MSE_GAMP_M, 3)./mean(Gnorm, 3)), 'r');
xlabel('SNR', 'fontsize',14)
ylabel('Normalized MSE', 'fontsize',14)
title('GAMP algorithm')

10*log10(mean(MSE_GAMP_M, 3)./mean(Gnorm, 3))

figure,
subplot(121)
gtrue2 = reshape(permute(reshape(gtrue,Nr,Nt,Nd),[1,3,2]),Nd*Nr,Nt);
bar3(abs(gtrue2));
title('True value of the channel', 'Fontsize', 14)
xlabel('Tx angle', 'Fontsize', 12)
ylabel('Rx angle & delay', 'Fontsize', 12)

subplot(122)
gtrue2 = reshape(permute(reshape(ghat_GAMP,Nr,Nt,Nd),[1,3,2]),Nd*Nr,Nt);
bar3(abs(gtrue2));
title('Estimate of the channel', 'Fontsize', 14)
xlabel('Tx angle', 'Fontsize', 12)
ylabel('Rx angle & delay', 'Fontsize', 12)



%%%% LS

figure,
plot(SNR_vector, 10*log10(mean(MSE_LS, 3)./mean(Gnorm, 3)), 'r');
xlabel('SNR', 'fontsize',14)
ylabel('Normalized MSE', 'fontsize',14)
title('Least Squares')

10*log10(mean(MSE_LS, 3)./mean(Gnorm, 3))

figure,
subplot(121)
gtrue2 = reshape(permute(reshape(gtrue,Nr,Nt,Nd),[1,3,2]),Nd*Nr,Nt);
bar3(abs(gtrue2));
title('True value of the channel', 'Fontsize', 14)
xlabel('Tx angle', 'Fontsize', 12)
ylabel('Rx angle & delay', 'Fontsize', 12)

subplot(122)
gtrue2 = reshape(permute(reshape(ghat_LS,Nr,Nt,Nd),[1,3,2]),Nd*Nr,Nt);
bar3(abs(gtrue2));
title('Estimate of the channel', 'Fontsize', 14)
xlabel('Tx angle', 'Fontsize', 12)
ylabel('Rx angle & delay', 'Fontsize', 12)


%% Delete the bad value
ii = 1;
for index_channel = 1:1:20
    if isequal (sort( reshape(MSE_LS(:,:, index_channel), length(SNR_vector), []), 'descend' ) ...
            , reshape(MSE_LS(:,:, index_channel), length(SNR_vector), []))
        MSE_LS_M(:,:, ii) = MSE_LS(:,:, index_channel);
        ii = ii + 1;
    end;
end;

Outlier_prob_LS = (ii-1)/20

figure,
plot(SNR_vector, 10*log10(mean(MSE_LS_M, 3)./mean(Gnorm, 3)), 'r');
xlabel('SNR', 'fontsize',14)
ylabel('Normalized MSE', 'fontsize',14)
title('GAMP algorithm')

10*log10(mean(MSE_GAMP, 3)./mean(Gnorm, 3))

figure,
subplot(121)
gtrue2 = reshape(permute(reshape(gtrue,Nr,Nt,Nd),[1,3,2]),Nd*Nr,Nt);
bar3(abs(gtrue2));
title('True value of the channel', 'Fontsize', 14)
xlabel('Tx angle', 'Fontsize', 12)
ylabel('Rx angle & delay', 'Fontsize', 12)

subplot(122)
gtrue2 = reshape(permute(reshape(ghat_LS,Nr,Nt,Nd),[1,3,2]),Nd*Nr,Nt);
bar3(abs(gtrue2));
title('Estimate of the channel', 'Fontsize', 14)
xlabel('Tx angle', 'Fontsize', 12)
ylabel('Rx angle & delay', 'Fontsize', 12)