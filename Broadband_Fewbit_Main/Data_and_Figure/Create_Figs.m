% clc;
clearvars;
dbstop if error;
close all;
rng(42)

fig_flag = 7;

Setup.N_sc = 16;
Setup.Nc = 512*20;
switch fig_flag
    case 1
        Setup.bit_vector = [1 2 3 4 5 6 7 8 +inf];
        Setup.Np_vector = 512*(2:1:10); %[251 509 1021 1531 2039 2557 3067 3583 4093];
        Setup.SNRdB_vector = linspace(-10, 40,11); %linspace(-10, 40, 11);
        Setup.Num_chan = 10;
        Setup.compute_rate_flag = 1;

        Setup.EM_BG_VAMP_flag = 1;
        Setup.EM_GM_VAMP_flag = 1;
        Setup.EM_BG_GAMP_flag = 1;
        Setup.EM_GM_GAMP_flag = 1;
        
        Setup.LS_flag = 1;
        Setup.LMMSE_flag = 1;
        Setup.BPDN_flag = 1;
        Setup.QIHT_flag = 1;
        
        [Results] = func_Broadband_Few_Bit_Comparison(Setup);
        save results_Chu_02252017 Setup Results
    case 2
        Setup.bit_vector = [1 2 3 4 5 6 7 8 +inf];
        Setup.Np_vector = 512*(2:1:10); %[251 509 1021 1531 2039 2557 3067 3583 4093];
        Setup.SNRdB_vector = linspace(-10, 40,11); %linspace(-10, 40, 11);
        Setup.Num_chan = 10;
        Setup.compute_rate_flag = 1;

        Setup.EM_BG_VAMP_flag = 1;
        Setup.EM_GM_VAMP_flag = 1;
        Setup.EM_BG_GAMP_flag = 1;
        Setup.EM_GM_GAMP_flag = 1;
        
        Setup.LS_flag = 1;
        Setup.LMMSE_flag = 1;
        Setup.BPDN_flag = 1;
        Setup.QIHT_flag = 1;
        
        [Results] = func_Broadband_Few_Bit_Comparison(Setup);
        save results_Chu_02252017 Setup Results      
        
    case 6
        Setup.bit_vector = 4;
        Setup.Np_vector = 2048; %[251 509 1021 1531 2039 2557 3067 3583 4093];
        Setup.SNRdB_vector = 0; %linspace(-10, 40,11); %linspace(-10, 40, 11);
        Setup.Num_chan = 1;
        Setup.compute_rate_flag = 0;

        Setup.EM_BG_VAMP_flag = 0;
        Setup.EM_GM_VAMP_flag = 1;
        Setup.EM_BG_GAMP_flag = 0;
        Setup.EM_GM_GAMP_flag = 0;
        
        Setup.LS_flag = 0;
        Setup.LMMSE_flag = 0;
        Setup.BPDN_flag = 0;
        Setup.QIHT_flag = 0;
        
        [Results] = func_Broadband_Few_Bit_Comparison(Setup);
%         save results_Chu_02252017 Setup Results  

case 7
        Setup.bit_vector = 4;
        Setup.Np_vector = 2048; %[251 509 1021 1531 2039 2557 3067 3583 4093];
        Setup.SNRdB_vector = 0; %linspace(-10, 40,11); %linspace(-10, 40, 11);
        Setup.Num_chan = 1;
        Setup.compute_rate_flag = 0;

        Setup.EM_BG_VAMP_flag = 0;
        Setup.EM_GM_VAMP_flag = 1;
        Setup.EM_BG_GAMP_flag = 0;
        Setup.EM_GM_GAMP_flag = 0;
        
        Setup.LS_flag = 0;
        Setup.LMMSE_flag = 0;
        Setup.BPDN_flag = 0;
        Setup.QIHT_flag = 0;
        
        [Results] = func_Broadband_Few_Bit_Comparison(Setup);
%         save results_Chu_02252017 Setup Results 
        
    case 8
        Setup.bit_vector = [1 2 3 4 +inf];
        Setup.Np_vector = 2048; %512*(2:1:10); %[251 509 1021 1531 2039 2557 3067 3583 4093];
        Setup.SNRdB_vector = linspace(-10, 40,11); %linspace(-10, 40, 11);
        Setup.Num_chan = 2;
        Setup.compute_rate_flag = 1;

        Setup.EM_BG_VAMP_flag = 1;
        Setup.EM_GM_VAMP_flag = 1;
        Setup.EM_BG_GAMP_flag = 1;
        Setup.EM_GM_GAMP_flag = 1;
        
        Setup.LS_flag = 0;
        Setup.LMMSE_flag = 0;
        Setup.BPDN_flag = 0;
        Setup.QIHT_flag = 0;
        
        [Results] = func_Broadband_Few_Bit_Comparison(Setup);
        
        save results_GAMP_VAMP_Chu_Np_2048_02252017.mat Results Setup
    
    case 9
        Setup.bit_vector = [1 2 3 4 5 6 7 8 +inf];
        Setup.Np_vector = 2048; %[251 509 1021 1531 2039 2557 3067 3583 4093];
        Setup.SNRdB_vector = 10; %linspace(-10, 40,11); %linspace(-10, 40, 11);
        Setup.Num_chan = 5;
        Setup.compute_rate_flag = 1;

        Setup.EM_BG_VAMP_flag = 1;
        Setup.EM_GM_VAMP_flag = 1;
        Setup.EM_BG_GAMP_flag = 1;
        Setup.EM_GM_GAMP_flag = 1;
        
        Setup.LS_flag = 1;
        Setup.LMMSE_flag = 1;
        Setup.BPDN_flag = 1;
        Setup.QIHT_flag = 1;
        
        [Results] = func_Broadband_Few_Bit_Comparison(Setup);
        save results_Chu_10db_Np_2048_bit_1_inf_02252017.mat Setup Results  
        
    case 13
        Setup.bit_vector = [1 2 3 4 +inf];
        Setup.Np_vector = 512*(2:1:10); %[251 509 1021 1531 2039 2557 3067 3583 4093];
        Setup.SNRdB_vector = 10; %linspace(-10, 40,11); %linspace(-10, 40, 11);
        Setup.Num_chan = 1;
        Setup.compute_rate_flag = 1;

        Setup.EM_BG_VAMP_flag = 1;
        Setup.EM_GM_VAMP_flag = 1;
        Setup.EM_BG_GAMP_flag = 1;
        Setup.EM_GM_GAMP_flag = 1;
        
        Setup.LS_flag = 0;
        Setup.LMMSE_flag = 0;
        Setup.BPDN_flag = 0;
        Setup.QIHT_flag = 0;
        
        [Results] = func_Broadband_Few_Bit_Comparison(Setup);
        
        save results_Chu_10dB_VAMP_GAMP_Np_1024_5120_02252017.mat Results Setup
        
    case 14
        Setup.bit_vector = 2;
        Setup.Np_vector = 512*(2:1:10); %[251 509 1021 1531 2039 2557 3067 3583 4093];
        Setup.SNRdB_vector = 10; %linspace(-10, 40, 6); %linspace(-10, 40, 11);
        Setup.Num_chan = 10;
        Setup.compute_rate_flag = 1;

        Setup.EM_BG_VAMP_flag = 1;
        Setup.EM_GM_VAMP_flag = 1;
        Setup.EM_BG_GAMP_flag = 1;
        Setup.EM_GM_GAMP_flag = 1;
        
        Setup.LS_flag = 1;
        Setup.LMMSE_flag = 1;
        Setup.BPDN_flag = 1;
        Setup.QIHT_flag = 1;
        [Results] = func_Broadband_Few_Bit_Comparison(Setup);
        save results_Chu_Rate_Np_2bit_10dB_02252017.mat Results Setup
end;