function [ Results ] = func_Broadband_Few_Bit_Comparison( Input )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% clc;
% clearvars;
% dbstop if error;
% close all;
% rng(42)

%% Model parameters
%% Channel Model: small test case
% load H_UPA_16_4_16_2clusters.mat
% Nta = 4;  % azimuth departure angle
% Nte = 4;  % elevation departure angle
% Nt = Nta * Nte;
% Nra = 2;  % azimuth arrival angle
% Nre = 2;  % elevation
% Nr = Nra * Nre;
% Nd = 16;  % Delay spread L
% debias = 1;

%% Channel Model: simulation setup
load H_UPA_64_64_16_4clusters.mat
Nta = 8;  % azimuth
Nte = 8;  % elevation
Nt = Nta * Nte;
Nra = 8;  % azimuth
Nre = 8;  % elevation
Nr = Nra * Nre;
Nd = 16;  % Delay spread L
debias = 1; %true;

%% Basis matrices
Br = 1/sqrt(Nr) * kron( dftmtx(Nra), dftmtx(Nre) );
Bt = 1/sqrt(Nt) * kron( dftmtx(Nta), dftmtx(Nte) );

%% ADC setup
bit_vector = Input.bit_vector; %2; %[1 2 3 4 +inf]; % [1, 2, 3, 4, 5, 6, 7, 8, +inf]; %[1 2 3 4 5 6 7 8 +inf];

% distortion facotor for the Additive Quantization Noise Model
% beta_vector = [0.3634 0.1175 0.03454 0.009497 0.002499 0.0006642 0.000166 0.00004151];
beta_vector = [0.3634 0.1188 0.03744 0.01154 0.003504 0.001035 0.0002999 0.00008543];

% SQNR = 10*log10(1./beta_vector);
% SQNR = [4.3962    9.2518   14.2666   19.3779   24.5544   29.8506
% 35.2302   40.6839];

%% stepsize
Lloyd_stepsize_vector = [1.5956 0.9957 0.586 0.3352 0.1881 0.1041 0.0569 0.0308];
stepsize_scale_factor = 1;

%% Pilot block length
Np_vector = Input.Np_vector; %2048; %512*(2:1:10); %[251 509 1021 1531 2039 2557 3067 3583 4093];

%% Training length
Nc = Input.Nc; %512*20;

%% SNR
SNRdB_vector = Input.SNRdB_vector; %linspace(-10, 40, 6); %linspace(-10, 40, 11);

%% number of channels for averaging
Num_chan = Input.Num_chan; %1;

%% number of subcarrier
N_sc = Input.N_sc; %16;

%% Choice of training sequences
% Sequences = 'Random_Gaussian';
% Sequences = 'Random_QPSK';
% Sequences = 'Shifted_Golay';
Sequences = 'Shifted_ZC';

% % % Sequences = 'Shifted_Pulse';
% % % Sequences = 'Shifted_QPSK';
% % % Sequences = 'Kronecker';
% Sequences = 'Golay';

Params.Sequences = Sequences;

if ( strcmp(Params.Sequences, 'Shifted_ZC')||...
        strcmp(Params.Sequences, 'Shifted_QPSK')||strcmp(Params.Sequences,'Shifted_Pulse')...
        || strcmp(Params.Sequences, 'Shifted_Golay'))
    Params.Enable_fast_implementation = 1; % Enable fast implementaion of Ax and A*z;
else
    Params.Enable_fast_implementation = 0; % Disable fast implementaion of Ax and A*z;
end;

%% struct Params contains the parameters of the system setup
Params.Nd = Nd;
Params.Nta = Nta;
Params.Nte = Nte;
Params.Nt = Nt;
Params.Nra = Nra;
Params.Nre = Nre;
Params.Nr = Nr;
Params.Nc = Nc;
Params.N_sc = N_sc;
Params.Br = Br;
Params.Bt = Bt;

Params.bit_vector = bit_vector;
Params.SNRdB_vector = SNRdB_vector;
Params.Np_vector = Np_vector;
Params.Num_chan = Num_chan;

%% run(or not run) algorithms
compute_rate_flag = Input.compute_rate_flag;

EM_BG_VAMP_flag = Input.EM_BG_VAMP_flag;
EM_GM_VAMP_flag = Input.EM_GM_VAMP_flag;
EM_BG_GAMP_flag = Input.EM_BG_GAMP_flag;
EM_GM_GAMP_flag = Input.EM_GM_GAMP_flag;

LS_flag = Input.LS_flag;
LMMSE_flag = Input.LMMSE_flag;
BPDN_flag = Input.BPDN_flag;
QIHT_flag = Input.QIHT_flag;

Lasso_flag = 0;
QIST_flag = 0;
CVTOSHT_flag = 0;

%% define the variables storing simulation results
MSE_BG_VAMP = zeros(length(SNRdB_vector), length(bit_vector), length(Np_vector), Num_chan);
MSE_GM_VAMP = zeros(length(SNRdB_vector), length(bit_vector), length(Np_vector), Num_chan);
MSE_BG_GAMP = zeros(length(SNRdB_vector), length(bit_vector), length(Np_vector), Num_chan);
MSE_GM_GAMP = zeros(length(SNRdB_vector), length(bit_vector), length(Np_vector), Num_chan);
MSE_LS      = zeros(length(SNRdB_vector), length(bit_vector), length(Np_vector), Num_chan);
MSE_LMMSE   = zeros(length(SNRdB_vector), length(bit_vector), length(Np_vector), Num_chan);
MSE_BPDN    = zeros(length(SNRdB_vector), length(bit_vector), length(Np_vector), Num_chan);
MSE_Lasso   = zeros(length(SNRdB_vector), length(bit_vector), length(Np_vector), Num_chan);
MSE_QIHT    = zeros(length(SNRdB_vector), length(bit_vector), length(Np_vector), Num_chan);
% MSE_QIST    = zeros(length(SNRdB_vector), length(bit_vector), length(Np_vector), Num_chan);
% MSE_CVTOSHT = zeros(length(SNRdB_vector), length(bit_vector), length(Np_vector), Num_chan);

NMSE_BG_VAMP_dB = zeros(length(SNRdB_vector), length(bit_vector), length(Np_vector));
NMSE_GM_VAMP_dB = zeros(length(SNRdB_vector), length(bit_vector), length(Np_vector));
NMSE_BG_GAMP_dB = zeros(length(SNRdB_vector), length(bit_vector), length(Np_vector));
NMSE_GM_GAMP_dB = zeros(length(SNRdB_vector), length(bit_vector), length(Np_vector));
NMSE_LS_dB      = zeros(length(SNRdB_vector), length(bit_vector), length(Np_vector));
NMSE_LMMSE_dB   = zeros(length(SNRdB_vector), length(bit_vector), length(Np_vector));
NMSE_BPDN_dB    = zeros(length(SNRdB_vector), length(bit_vector), length(Np_vector));
NMSE_Lasso_dB   = zeros(length(SNRdB_vector), length(bit_vector), length(Np_vector));
NMSE_QIHT_dB    = zeros(length(SNRdB_vector), length(bit_vector), length(Np_vector));
% NMSE_QIST_dB    = zeros(length(SNRdB_vector), length(bit_vector), length(Np_vector));
% NMSE_CVTOSHT_dB = zeros(length(SNRdB_vector), length(bit_vector), length(Np_vector));

Rate_pcsi   = zeros(length(SNRdB_vector), length(bit_vector), length(Np_vector), Num_chan);
Rate_BG_VAMP= zeros(length(SNRdB_vector), length(bit_vector), length(Np_vector), Num_chan);
Rate_GM_VAMP= zeros(length(SNRdB_vector), length(bit_vector), length(Np_vector), Num_chan);
Rate_BG_GAMP= zeros(length(SNRdB_vector), length(bit_vector), length(Np_vector), Num_chan);
Rate_GM_GAMP= zeros(length(SNRdB_vector), length(bit_vector), length(Np_vector), Num_chan);
Rate_LS     = zeros(length(SNRdB_vector), length(bit_vector), length(Np_vector), Num_chan);
Rate_LMMSE  = zeros(length(SNRdB_vector), length(bit_vector), length(Np_vector), Num_chan);
Rate_BPDN   = zeros(length(SNRdB_vector), length(bit_vector), length(Np_vector), Num_chan);
Rate_Lasso  = zeros(length(SNRdB_vector), length(bit_vector), length(Np_vector), Num_chan);
Rate_QIHT   = zeros(length(SNRdB_vector), length(bit_vector), length(Np_vector), Num_chan);
% Rate_QIST   = zeros(length(SNRdB_vector), length(bit_vector), length(Np_vector), Num_chan);
% Rate_CVTOSHT= zeros(length(SNRdB_vector), length(bit_vector), length(Np_vector), Num_chan);

% Rate_pcsi_mean   = zeros(length(SNRdB_vector), length(bit_vector), length(Np_vector));
Rate_BG_VAMP_mean= zeros(length(SNRdB_vector), length(bit_vector), length(Np_vector));
Rate_GM_VAMP_mean= zeros(length(SNRdB_vector), length(bit_vector), length(Np_vector));
Rate_BG_GAMP_mean= zeros(length(SNRdB_vector), length(bit_vector), length(Np_vector));
Rate_GM_GAMP_mean= zeros(length(SNRdB_vector), length(bit_vector), length(Np_vector));
Rate_LS_mean     = zeros(length(SNRdB_vector), length(bit_vector), length(Np_vector));
Rate_LMMSE_mean  = zeros(length(SNRdB_vector), length(bit_vector), length(Np_vector));
Rate_BPDN_mean   = zeros(length(SNRdB_vector), length(bit_vector), length(Np_vector));
Rate_QIHT_mean   = zeros(length(SNRdB_vector), length(bit_vector), length(Np_vector));
% Rate_Lasso_mean  = zeros(length(SNRdB_vector), length(bit_vector), length(Np_vector));
% Rate_QIST_mean   = zeros(length(SNRdB_vector), length(bit_vector), length(Np_vector));
% Rate_CVTOSHT_mean= zeros(length(SNRdB_vector), length(bit_vector), length(Np_vector));

MI_pcsi     = zeros(length(SNRdB_vector), length(bit_vector), length(Np_vector), Num_chan);
MI_BG_VAMP  = zeros(length(SNRdB_vector), length(bit_vector), length(Np_vector), Num_chan);
MI_GM_VAMP  = zeros(length(SNRdB_vector), length(bit_vector), length(Np_vector), Num_chan);
MI_BG_GAMP  = zeros(length(SNRdB_vector), length(bit_vector), length(Np_vector), Num_chan);
MI_GM_GAMP  = zeros(length(SNRdB_vector), length(bit_vector), length(Np_vector), Num_chan);
MI_LS       = zeros(length(SNRdB_vector), length(bit_vector), length(Np_vector), Num_chan);
MI_LMMSE    = zeros(length(SNRdB_vector), length(bit_vector), length(Np_vector), Num_chan);
MI_BPDN     = zeros(length(SNRdB_vector), length(bit_vector), length(Np_vector), Num_chan);
MI_QIHT     = zeros(length(SNRdB_vector), length(bit_vector), length(Np_vector), Num_chan);
% MI_Lasso    = zeros(length(SNRdB_vector), length(bit_vector), length(Np_vector), Num_chan);
% MI_QIST     = zeros(length(SNRdB_vector), length(bit_vector), length(Np_vector), Num_chan);
% MI_CVTOSHT  = zeros(length(SNRdB_vector), length(bit_vector), length(Np_vector), Num_chan);

MI_pcsi_mean     = zeros(length(SNRdB_vector), length(bit_vector), length(Np_vector));
MI_BG_VAMP_mean  = zeros(length(SNRdB_vector), length(bit_vector), length(Np_vector));
MI_GM_VAMP_mean  = zeros(length(SNRdB_vector), length(bit_vector), length(Np_vector));
MI_BG_GAMP_mean  = zeros(length(SNRdB_vector), length(bit_vector), length(Np_vector));
MI_GM_GAMP_mean  = zeros(length(SNRdB_vector), length(bit_vector), length(Np_vector));
MI_LS_mean       = zeros(length(SNRdB_vector), length(bit_vector), length(Np_vector));
MI_LMMSE_mean    = zeros(length(SNRdB_vector), length(bit_vector), length(Np_vector));
MI_BPDN_mean     = zeros(length(SNRdB_vector), length(bit_vector), length(Np_vector));
MI_QIHT_mean     = zeros(length(SNRdB_vector), length(bit_vector), length(Np_vector));
% MI_Lasso_mean    = zeros(length(SNRdB_vector), length(bit_vector), length(Np_vector), Num_chan);
% MI_QIST     = zeros(length(SNRdB_vector), length(bit_vector), length(Np_vector), Num_chan);
% MI_CVTOSHT  = zeros(length(SNRdB_vector), length(bit_vector), length(Np_vector), Num_chan);

time_BG_VAMP_vector = [];
time_GM_VAMP_vector = [];
time_BG_GAMP_vector = [];
time_GM_GAMP_vector = [];
time_BPDN_vector = [];

Gnorm = nan(1, Num_chan);

for index_channel = 1:1:Num_chan
    % virtual channel model or angular domain channel model
    %         Gtrue = func_mmWave_Channel_UPA(Nt_W, Nt_H, Nr_W, Nr_H, ...
    %             Angle_spread_W, Angle_spread_H, Num_cluster, Num_ray, Narrowband_flag, Nd);
    Gtrue = Hv_all(:,:,:,index_channel);
    %     Gtrue = Hv_all(:,:,:,3);
    gtrue = reshape(Gtrue, [],1);
    
    Gnorm(index_channel) = norm(gtrue);
    
    for index_Np = 1:1:length(Np_vector)
        
        Np = Np_vector(index_Np);
        Params.Np = Np;
        
        % Training sequences in the spatial and time domain, size: Nt X Np
        % Walsh_Hadamard_matrix = hadamard(Np);
        % T_train = Walsh_Hadamard_matrix(randsample(1:Np,Nt), 1:Np);
        
        switch Params.Sequences
            case 'Random_QPSK'                 % random QPSK sequence
                T_train_unit = sign(randn(Nt,Np))/sqrt(2)+1i*sign(randn(Nt,Np))/sqrt(2);
                
            case 'Random_Gaussian'
                T_train_unit = randn(Nt,Np)/sqrt(2) + 1i*randn(Nt,Np)/sqrt(2);
                
            case 'Random_ZC'
                T_train_unit = func_ZadoffChuSeq(Nt, Np);
                
            case 'Golay'
                if floor(log2(Np)) == log2(Np)
                    [a, b] = generate_golay(log2(Np)-1);
                else
                    keyboard;
                end;
                s = [a, b];
                T0 = zeros(Nt,Np);
                for i=0:Nt-1
                    %                     T0(i+1,:) = circshift(s,[0 i*Nd]);
                    T0(i+1,:) = [a(Np/2-i*Nd/2+1: Np/2), a(1:Np/2-i*Nd/2), ...
                        b(Np/2-i*Nd/2+1: Np/2), b(1:Np/2-i*Nd/2)];
                    %for this kind of T0, the singular values are same
                end
                T_train_unit = T0;
                
            case 'Shifted_ZC'
                % This training sequence design is proposed by Prof. Schniter
                if mod(Np,2) == 0
                    s = exp((1i*pi/Np)*(0:Np-1).^2);  % even length
                else
                    s = exp((1i*pi/Np)*(0:Np-1).*((0:Np-1)+1)); % odd length
                end;
                T0 = zeros(Nt,Np);
                for i=0:Nt-1
                    T0(i+1,:) = circshift(s,[0 i*Nd]);
                end
                T_train_unit = T0;
                
            case 'Shifted_QPSK'
                s = sign(randn(1,Np))/sqrt(2)+1i*sign(randn(1,Np))/sqrt(2);
                T0 = zeros(Nt,Np);
                for i=0:Nt-1
                    T0(i+1,:) = circshift(s,[0 i*Nd]);
                end
                T_train_unit = T0;
                
            case 'Kronecker'
                % [1 0 0 0 ...]
                % This training sequence design is proposed by Prof. Schniter
                s = zeros(1, Np);
                s(1,1) = 1;
                T0 = zeros(Nt,Np);
                for i=0:Nt-1
                    T0(i+1,:) = circshift(s,[0 i*Nd]);
                end
                T_train_unit = T0;
                
            case  'Shifted_Pulse'
                s = zeros(Nd, Np/Nd);
                if mod(Np/Nd,2) == 0
                    s(1,:) = exp((1i*pi/ (Np/Nd) )*(0:Np/Nd-1).^2);  % even length
                else
                    s(1,:) = exp((1i*pi/ (Np/Nd) )*(0:Np/Nd-1).*((0:Np/Nd-1)+1)); % odd length
                end;
                s = reshape(s, 1, Np);
                T0 = zeros(Nt,Np);
                for i=0:Nt-1
                    T0(i+1,:) = circshift(s,[0 i*Nd]);
                end
                T_train_unit = T0;
                
            case  'Shifted_Golay'
                display('Not well defined!');
                keyboard;
                %                 if floor(log2(Np)) == log2(Np)
                %                     [a, b] = generate_golay(log2(Np)-1);
                %                 else
                %                     keyboard;
                %                 end;
                %                 s = [a, b];
                %                 T0 = zeros(Nt,Np);
                %                 for i=0:Nt-1
                %                     T0(i+1,:) = circshift(s,[0 i*Nd]);
                %                 end
                %                 T_train_unit = T0;
        end
        
        
        %                 T_cyclic_0 = NaN(Nr*Nd, Np);
        %                 for ii=0:1:Nd-1
        %                     T_cyclic_0( ii*Nt+1: (ii+1)*Nt, :) = [T_train_scaled(:,Np-ii+1:Np), T_train_scaled(:,1:Np-ii)];
        %                 end
        
        %% fast construction of T_cyclic by matlab function 'circshift'
        T_cyclic_unit = NaN(Nr*Nd, Np);
        for ii=0:1:Nd-1
            T_cyclic_unit( ii*Nt+1: (ii+1)*Nt, :) = circshift(T_train_unit,ii,2);
        end
        
        % % 1st fast implementation of V_cyclic = kron(eye(Nd), Bt') * T_cyclic
        %                 T_temp = reshape(T_cyclic, Nt, Nd, Np);
        %                 V_temp = Bhfast(T_temp, Nt);
        %                 V_cyclic_1 = reshape(V_temp, Nd*Nt, Np);
        
        % % 2nd fast implementation of V_cyclic = kron(eye(Nd),
        % Bt') * T_cyclic. The time consumptions of 1st and 2nd
        % implementation are similar.
        V_cyclic_unit = NaN(Nr*Nd, Np);
        V_train_unit = Bhfast(T_train_unit, Nt);
        for ii=0:1:Nd-1
            V_cyclic_unit( ii*Nt+1: (ii+1)*Nt, :) = circshift(V_train_unit,ii,2);
        end
        
        % compute V (right singular matrix) of the measurement matrix A
        if ~strcmp(Params.Sequences, 'Shifted_ZC')
            VcVt_unit = conj(V_cyclic_unit) * transpose(V_cyclic_unit);
            [V_of_VcVt_unit, D_of_VcVt_unit] = eig(0.5*(VcVt_unit + ctranspose(VcVt_unit)));
        end;
        
        for index_bit = 1:1:length(bit_vector)
            
            bit = bit_vector(index_bit);
            if bit~=+inf
                beta = beta_vector(bit);
            else
                beta = 0;
            end;
            
            display(index_channel);
            display(Np);
            display(bit);
            for index_SNR = 1:1:length(SNRdB_vector)
                
                SNRdB = SNRdB_vector(index_SNR);
                
                SNR_linear = 10^(SNRdB/10);
                
                SNR_coeff = 1/(norm(T_train_unit, 'fro'))* sqrt(Np) * sqrt(SNR_linear);
                
                T_cyclic  = T_cyclic_unit.* SNR_coeff;
                % scale T_train such that sum((T_train_scaled(:,ii)).^2) = SNR_linear for
                % ii=1,2,...
                V_cyclic = V_cyclic_unit.* SNR_coeff;
                Params.V_cyclic = V_cyclic;
                if ~strcmp(Params.Sequences, 'Shifted_ZC')
                    Params.V_of_VcVt =  V_of_VcVt_unit;
                    Params.D_of_VcVt = diag(D_of_VcVt_unit).*SNR_coeff^2;
                end;
                
                % z = A *gtrue;      % "True" transform vector
                %                 z_0 = reshape( Br * reshape(Gtrue, Nr, []) * kron(eye(Nd), Bt') * T_cyclic, [], 1);
                z = reshape( Br * reshape(Gtrue, Nr, []) * V_cyclic, [], 1);
                
                wvar = 1; % the noise variance is fixed to be one
                
                y = z + sqrt(1/2)*wvar*randn(Nr*Np,2)*[1;1i];
                
                %% norm esimation (seems to be the MMSE estimator??)
                % compute the variance gain due to the measurement matrix 'A'
                %                 gain = trace(T_cyclic * T_cyclic')/Np;
                gain = SNR_linear * Nd;
                
                % estimate the variance of each element of gtrue
                gvar = ( mean(abs(y).^2) - wvar )/gain;  % is close to SNR_linear/gain
                Params.gvar = gvar; % is close to 1/Nd
                
                % estimate the norm of gtrue (note the dimension of gtrue is Nt*Nr*Nd)
                estimated_norm = sqrt(gvar * Nt * Nr * Nd);
                Params.estimated_norm = estimated_norm;
                %                     display(estimated_norm);
                %                     gture_norm = norm(gtrue);
                
                %% quantization
                
                %% dithering
                dither_mean = rms([real(y); imag(y)])*0.0; % default value of dithering threshold
                
                %% ADC resolution
                if bit == +inf
                    r_bit = y;
                    stepsize = 0;
                    beta = 0;
                    %                     r = r_bit;
                else
                    stepsize = stepsize_scale_factor * rms([real(y); imag(y)])* Lloyd_stepsize_vector(bit);
                    %                         stepsize = ( 1.01 * max(abs([real(y); imag(y)])))/(2^(bit-1)); % previous code
                    
                    r_bit_real = min ( max( floor( (real(y)-dither_mean)/stepsize) + 2^(bit - 1), 0), 2^bit-1) ; %[0 2^bit-1]
                    r_bit_imag = min ( max( floor( (imag(y)-dither_mean)/stepsize )+ 2^(bit - 1), 0), 2^bit-1) ; %[0 2^bit-1]
                    
                    r_bit = r_bit_real + 1j * r_bit_imag;
                    
                    %                     r = sign(real(y)) .* ( min( ceil( abs(real(y)) /stepsize) , 2^(bit-1) ) - 1/2 ) * stepsize  + ...
                    %                         1j* sign(imag(y)) .* ( min( ceil( abs(imag(y)) /stepsize) , 2^(bit-1) ) - 1/2 ) * stepsize;
                end;
                
                % for the testing of BPDN
                %                 mu = mean(real(y).*real(r))/rms(real(y))^2
                %                 sigma = sqrt ( rms(( real(r) - mu.* real(y)))^2 + mu^2 * 1/2 )
                
                Params.bit = bit;
                Params.stepsize = stepsize;
                Params.beta = beta;
                Params.SNR_linear = SNR_linear;
                Params.dither_mean = dither_mean;
                Params.wvar = wvar;
                
                % only used in the BPDN algorithm
                Params.stepsize_scale_factor = stepsize_scale_factor;
                Params.received_power_linear = rms([real(y); imag(y)])^2;
                
                % only used for debugging
                Params.Gtrue = Gtrue;
                
                % EM-BG-VAMP algorithm
                if EM_BG_VAMP_flag == 1
                    rng(index_channel);
                    
                    [ ghat_BG_VAMP, time_VAMP]   = func_Broadband_Few_Bit_UPA_EMBGVAMP( T_cyclic, r_bit, Params);
                    time_BG_VAMP_vector = [time_BG_VAMP_vector, time_VAMP];
                    Ghat_BG_VAMP = reshape(ghat_BG_VAMP, Nr, Nt, Nd);
                    if debias
                        ghat_BG_VAMP = estimated_norm * ghat_BG_VAMP /norm(ghat_BG_VAMP);
                    end;
                    MSE_BG_VAMP(index_SNR, index_bit, index_Np, index_channel) = norm(gtrue - ghat_BG_VAMP, 2)^2;
                end;
                
                % EM-GM-VAMP algorithm
                if EM_GM_VAMP_flag == 1
                    rng(index_channel);
                    
                    [ ghat_GM_VAMP, time_VAMP]   = func_Broadband_Few_Bit_UPA_EMGMVAMP( T_cyclic, r_bit, Params);
                    time_GM_VAMP_vector = [time_GM_VAMP_vector, time_VAMP];
                    Ghat_GM_VAMP = reshape(ghat_GM_VAMP, Nr, Nt, Nd);
                    if debias
                        ghat_GM_VAMP = estimated_norm * ghat_GM_VAMP /norm(ghat_GM_VAMP);
                    end;
                    MSE_GM_VAMP(index_SNR, index_bit, index_Np, index_channel) = norm(gtrue - ghat_GM_VAMP, 2)^2;
                end;
                
                
                % EMBGAMP algorith
                if EM_BG_GAMP_flag == 1
                    rng(index_channel);
                    [ ghat_BG_GAMP,  time_BG]   = func_Broadband_Few_Bit_UPA_EMBGGAMP( T_cyclic, r_bit, Params);
                    time_BG_GAMP_vector = [time_BG_GAMP_vector, time_BG];
                    Ghat_BG_GAMP = reshape(ghat_BG_GAMP, Nr, Nt, Nd);
                    if debias
                        ghat_BG_GAMP = estimated_norm * ghat_BG_GAMP /norm(ghat_BG_GAMP);
                    end;
                    MSE_BG_GAMP(index_SNR, index_bit, index_Np, index_channel) = norm(gtrue - ghat_BG_GAMP, 2)^2;
                    %if debias, gain = (ghat_BG'*gtrue)/norm(ghat_BG)^2; else gain = 1; end;
                    %MSE_BG(index_SNR, index_bit, index_Np, index_channel) = norm(gtrue - gain*ghat_BG, 2)^2;
                    %                     GI_MSE_BG(index_SNR, index_bit, index_Np, index_channel) = norm( gtrue/norm(gtrue) - ghat_BG/norm(ghat_BG) )^2;
                end;
                
                % EMGMAMP algorithm
                if EM_GM_GAMP_flag == 1
                    rng(index_channel);
                    [ ghat_GM_GAMP, time_GM]   = func_Broadband_Few_Bit_UPA_EMGMGAMP( T_cyclic, r_bit, Params);
                    time_GM_GAMP_vector = [time_GM_GAMP_vector, time_GM];
                    Ghat_GM_GAMP = reshape(ghat_GM_GAMP, Nr, Nt, Nd);
                    if debias
                        ghat_GM_GAMP = estimated_norm * ghat_GM_GAMP /norm(ghat_GM_GAMP);
                    end;
                    MSE_GM_GAMP(index_SNR, index_bit, index_Np, index_channel) = norm(gtrue - ghat_GM_GAMP, 2)^2;
                    
                    %                     if debias, gain = (ghat_GM'*gtrue)/norm(ghat_GM)^2; else gain = 1; end;
                    %                     MSE_GM(index_SNR, index_bit, index_Np, index_channel) = norm(gtrue - gain*ghat_GM, 2)^2;
                    %                     GI_MSE_GM(index_SNR, index_bit, index_Np, index_channel) = norm( gtrue/norm(gtrue) - ghat_GM/norm(ghat_GM) )^2;
                end;
                
                
                % Least Squares algorithm
                if LS_flag == 1
                    rng(index_channel);
                    [ ghat_LS ] = func_Broadband_Few_Bit_UPA_LS( T_cyclic,  r_bit, Params);
                    Ghat_LS = reshape(ghat_LS, Nr, Nt, Nd);
                    if debias
                        ghat_LS = estimated_norm * ghat_LS /norm(ghat_LS);
                    end;
                    MSE_LS(index_SNR, index_bit, index_Np, index_channel) = norm(gtrue - ghat_LS, 2)^2;
                    %                     if debias, gain = (ghat_LS'*gtrue)/norm(ghat_LS)^2; else gain = 1; end;
                    %                     MSE_LS(index_SNR, index_bit, index_Np, index_channel) = norm(gtrue - gain*ghat_LS, 2)^2;
                    %                     GI_MSE_LS(index_SNR, index_bit, index_Np, index_channel) = norm( gtrue/norm(gtrue) - ghat_LS/norm(ghat_LS) )^2;
                end;
                
                % Linear MMSE algorithm
                if LMMSE_flag == 1
                    rng(index_channel);
                    [ ghat_LMMSE ] = func_Broadband_Few_Bit_UPA_ALMMSE( T_cyclic,  r_bit, Params);
                    Ghat_LMMSE = reshape(ghat_LMMSE, Nr, Nt, Nd);
                    if debias
                        ghat_LMMSE = estimated_norm * ghat_LMMSE /norm(ghat_LMMSE);
                    end;
                    MSE_LMMSE(index_SNR, index_bit, index_Np, index_channel) = norm(gtrue - ghat_LMMSE, 2)^2;
                    %                     if debias, gain = (ghat_LMMSE'*gtrue)/norm(ghat_LMMSE)^2; else gain = 1; end;
                    %                     MSE_LMMSE(index_SNR, index_bit, index_Np, index_channel) = norm(gtrue - gain*ghat_LMMSE, 2)^2;
                    %                     GI_MSE_LMMSE(index_SNR, index_bit, index_Np, index_channel) = norm( gtrue/norm(gtrue) - ghat_LMMSE/norm(ghat_LMMSE) )^2;
                end;
                
                % BPDN algorithm
                if BPDN_flag == 1
                    rng(index_channel);
                    [ ghat_BPDN,  time_BPDN ] = func_Broadband_Few_Bit_UPA_BPDN(T_cyclic, r_bit, Params);
                    time_BPDN_vector= [time_BPDN_vector, time_BPDN];
                    Ghat_BPDN = reshape(ghat_BPDN, Nr, Nt, Nd);
                    if debias
                        ghat_BPDN = estimated_norm * ghat_BPDN /norm(ghat_BPDN);
                    end;
                    MSE_BPDN(index_SNR, index_bit, index_Np, index_channel) = norm(gtrue - ghat_BPDN, 2)^2;
                    %                     if debias, gain = (ghat_BPDN'*gtrue)/norm(ghat_BPDN)^2; else gain = 1; end;
                    %                     MSE_BPDN(index_SNR, index_bit, index_Np, index_channel) = norm(gtrue - gain*ghat_BPDN, 2)^2;
                    %                     GI_MSE_BPDN(index_SNR, index_bit, index_Np, index_channel) = norm( gtrue/norm(gtrue) - ghat_BPDN/norm(ghat_BPDN) )^2;
                end;
                
                % Lasso algorithm
                if Lasso_flag == 1
                    rng(index_channel);
                    [ ghat_Lasso ] = func_Broadband_Few_Bit_UPA_Lasso(T_cyclic, r_bit, Params);
                    Ghat_Lasso = reshape(ghat_Lasso, Nr, Nt, Nd);
                    if debias
                        ghat_Lasso = estimated_norm * ghat_Lasso /norm(ghat_Lasso);
                    end;
                    MSE_Lasso(index_SNR, index_bit, index_Np, index_channel) = norm(gtrue - ghat_Lasso, 2)^2;
                    %                     if debias, gain = (ghat_BPDN'*gtrue)/norm(ghat_BPDN)^2; else gain = 1; end;
                    %                     MSE_BPDN(index_SNR, index_bit, index_Np, index_channel) = norm(gtrue - gain*ghat_BPDN, 2)^2;
                    %                     GI_MSE_BPDN(index_SNR, index_bit, index_Np, index_channel) = norm( gtrue/norm(gtrue) - ghat_BPDN/norm(ghat_BPDN) )^2;
                end;
                
                
                
                % QIST algorithm
                if QIST_flag == 1
                    rng(index_channel);
                    [ ghat_QIST ] = func_Broadband_Few_Bit_UPA_QIST(T_cyclic, r_bit, Params);
                    Ghat_QIST = reshape(ghat_QIST, Nr, Nt, Nd);
                    if debias
                        ghat_QIST = estimated_norm * ghat_QIST /norm(ghat_QIST);
                    end;
                    MSE_QIST(index_SNR, index_bit, index_Np, index_channel) = norm(gtrue - ghat_QIST, 2)^2;
                    %                     if debias, gain = (ghat_BPDN'*gtrue)/norm(ghat_BPDN)^2; else gain = 1; end;
                    %                     MSE_BPDN(index_SNR, index_bit, index_Np, index_channel) = norm(gtrue - gain*ghat_BPDN, 2)^2;
                    %                     GI_MSE_BPDN(index_SNR, index_bit, index_Np, index_channel) = norm( gtrue/norm(gtrue) - ghat_BPDN/norm(ghat_BPDN) )^2;
                end;
                
                % QIHT algorithm
                if QIHT_flag == 1
                    rng(index_channel);
                    [ ghat_QIHT ] = func_Broadband_Few_Bit_UPA_QIHT(T_cyclic, r_bit, Params);
                    Ghat_QIHT = reshape(ghat_QIHT, Nr, Nt, Nd);
                    if debias
                        ghat_QIHT = estimated_norm * ghat_QIHT /norm(ghat_QIHT);
                    end;
                    MSE_QIHT(index_SNR, index_bit, index_Np, index_channel) = norm(gtrue - ghat_QIHT, 2)^2;
                    %                     if debias, gain = (ghat_BPDN'*gtrue)/norm(ghat_BPDN)^2; else gain = 1; end;
                    %                     MSE_BPDN(index_SNR, index_bit, index_Np, index_channel) = norm(gtrue - gain*ghat_BPDN, 2)^2;
                    %                     GI_MSE_BPDN(index_SNR, index_bit, index_Np, index_channel) = norm( gtrue/norm(gtrue) - ghat_BPDN/norm(ghat_BPDN) )^2;
                end;
                
                % CVTOSHT algorithm
                if CVTOSHT_flag == 1
                    rng(index_channel);
                    [ ghat_CVTOSHT ] = func_Broadband_Few_Bit_UPA_CVTOSHT(T_cyclic,  r_bit, Params);
                    Ghat_CVTOSHT = reshape(ghat_CVTOSHT, Nr, Nt, Nd);
                    if debias
                        ghat_CVTOSHT = estimated_norm * ghat_CVTOSHT /norm(ghat_CVTOSHT);
                    end;
                    MSE_CVTOSHT(index_SNR, index_bit, index_Np, index_channel) = norm(gtrue - ghat_CVTOSHT, 2)^2;
                    %                     if debias, gain = (ghat_CVTOSHT'*gtrue)/norm(ghat_CVTOSHT)^2; else gain = 1; end;
                    %                     MSE_CVTOSHT(index_SNR, index_bit, index_Np, index_channel) = norm(gtrue - gain*ghat_CVTOSHT, 2)^2;
                    %                     GI_MSE_CVTOSHT(index_SNR, index_bit, index_Np, index_channel) = norm( gtrue/norm(gtrue) - ghat_CVTOSHT/norm(ghat_CVTOSHT) )^2;
                end;
                
                
                if compute_rate_flag == 1
                    % Rate computation with perfect CSI
                    % Should use SVD decomposition here
                    % the variance of noise
                    [MI, Rate] = eval_chanest(Gtrue,Gtrue, Params);
                    MI_pcsi(index_SNR, index_bit, index_Np, index_channel) = MI;
                    Rate_pcsi(index_SNR, index_bit, index_Np, index_channel) = Rate;
                    
                    if EM_BG_VAMP_flag == 1
                        [MI, Rate] = eval_chanest(Ghat_BG_VAMP,Gtrue, Params);
                        MI_BG_VAMP(index_SNR, index_bit, index_Np, index_channel) = MI;
                        Rate_BG_VAMP(index_SNR, index_bit, index_Np, index_channel) = Rate;
                    end;
                    
                    if EM_GM_VAMP_flag == 1
                        [MI, Rate] = eval_chanest(Ghat_GM_VAMP,Gtrue, Params);
                        MI_GM_VAMP(index_SNR, index_bit, index_Np, index_channel) = MI;
                        Rate_GM_VAMP(index_SNR, index_bit, index_Np, index_channel) = Rate;
                    end;
                    
                    if EM_BG_GAMP_flag == 1
                        [MI, Rate] = eval_chanest(Ghat_BG_GAMP,Gtrue, Params);
                        MI_BG_GAMP(index_SNR, index_bit, index_Np, index_channel) = MI;
                        Rate_BG_GAMP(index_SNR, index_bit, index_Np, index_channel) = Rate;
                    end;
                    
                    if EM_GM_GAMP_flag == 1
                        [MI, Rate] = eval_chanest(Ghat_GM_GAMP,Gtrue, Params);
                        MI_GM_GAMP(index_SNR, index_bit, index_Np, index_channel) = MI;
                        Rate_GM_GAMP(index_SNR, index_bit, index_Np, index_channel) = Rate;
                    end;
                    
                    % Rate computation with LS estimation of the CSI
                    if LS_flag == 1
                        [MI, Rate] = eval_chanest(Ghat_LS,Gtrue, Params);
                        MI_LS(index_SNR, index_bit, index_Np, index_channel) = MI;
                        Rate_LS(index_SNR, index_bit, index_Np, index_channel) = Rate;
                    end;
                    
                    
                    % Rate computation with LMMSE estimation of the CSI
                    if LMMSE_flag == 1
                        [MI, Rate] = eval_chanest(Ghat_LMMSE,Gtrue, Params);
                        MI_LMMSE(index_SNR, index_bit, index_Np, index_channel) = MI;
                        Rate_LMMSE(index_SNR, index_bit, index_Np, index_channel) = Rate;
                    end;
                    
                    % Rate computation with BPDN estimation of the CSI
                    if BPDN_flag == 1
                        [MI, Rate] = eval_chanest(Ghat_BPDN,Gtrue, Params);
                        MI_BPDN(index_SNR, index_bit, index_Np, index_channel) = MI;
                        Rate_BPDN(index_SNR, index_bit, index_Np, index_channel) = Rate;
                    end;
                    
                    % Rate computation with Lasso estimation of the CSI
                    if Lasso_flag == 1
                        [MI, Rate] = eval_chanest(Ghat_Lasso,Gtrue, Params);
                        MI_Lasso(index_SNR, index_bit, index_Np, index_channel) = MI;
                        Rate_Lasso(index_SNR, index_bit, index_Np, index_channel) = Rate;
                    end;
                    
                    % Rate computation with QIHT estimation of the CSI
                    if QIHT_flag == 1
                        [MI, Rate] = eval_chanest(Ghat_QIHT,Gtrue, Params);
                        MI_QIHT(index_SNR, index_bit, index_Np, index_channel) = MI;
                        Rate_QIHT(index_SNR, index_bit, index_Np, index_channel) = Rate;
                    end;
                    
                    % Rate computation with QIST estimation of the CSI
                    if QIST_flag == 1
                        [MI, Rate] = eval_chanest(Ghat_QIST,Gtrue, Params);
                        MI_QIST(index_SNR, index_bit, index_Np, index_channel) = MI;
                        Rate_QIST(index_SNR, index_bit, index_Np, index_channel) = Rate;
                    end;
                end;
            end;
        end;
    end;
end;


%% EM-BG-VAMP
if EM_BG_VAMP_flag == 1
    %%%% VAMP
    NMSE_BG_VAMP_dB = 10*log10(mean(MSE_BG_VAMP, 4)/mean(Gnorm.^2));
    % MSE_VAMP_dB = 10*log10(MSE_GAMP(:,:,:,4)./Gnorm(:,:,:,end));
    %     GI_MSE_BG_dB = 10*log10(mean(GI_MSE_BG, 4));
    
    %     % MSE versus SNR
    %     figure,
    %     plot(SNRdB_vector, NMSE_BG_VAMP_dB(:,:, end), '.-');
    %     xlabel('SNR', 'fontsize',14)
    %     ylabel('Normalized MSE', 'fontsize',14)
    %     if debias, ylabel('Normalized debiased MSE', 'fontsize',14); end;
    %     title('EM-BG-VAMP algorithm')
    %     grid on;
    %     legend('1-bit', '2-bit', '3-bit', '4-bit', '\infty-bit')
end;

%% EM-GM-VAMP
if EM_GM_VAMP_flag == 1
    %%%% VAMP
    NMSE_GM_VAMP_dB = 10*log10(mean(MSE_GM_VAMP, 4)/mean(Gnorm.^2));
    % MSE_VAMP_dB = 10*log10(MSE_GAMP(:,:,:,4)./Gnorm(:,:,:,end));
    %     GI_MSE_BG_dB = 10*log10(mean(GI_MSE_BG, 4));
    
    %     % MSE versus SNR
    %     figure,
    %     plot(SNRdB_vector, NMSE_GM_VAMP_dB(:,:, end), '.-');
    %     xlabel('SNR', 'fontsize',14)
    %     ylabel('Normalized MSE', 'fontsize',14)
    %     if debias, ylabel('Normalized debiased MSE', 'fontsize',14); end;
    %     title('EM-GM-VAMP algorithm')
    %     grid on;
    %     legend('1-bit', '2-bit', '3-bit', '4-bit', '\infty-bit')
end;

%% EM-BG-GAMP
if EM_BG_GAMP_flag == 1
    %%%% GAMP
    NMSE_BG_GAMP_dB = 10*log10(mean(MSE_BG_GAMP, 4)/mean(Gnorm.^2));
    % MSE_GAMP_dB = 10*log10(MSE_GAMP(:,:,:,4)./Gnorm(:,:,:,end));
    %     GI_MSE_BG_dB = 10*log10(mean(GI_MSE_BG, 4));
    
    %     % MSE versus SNR
    %     figure,
    %     plot(SNRdB_vector, NMSE_BG_GAMP_dB(:,:, end), '.-');
    %     xlabel('SNR', 'fontsize',14)
    %     ylabel('Normalized MSE', 'fontsize',14)
    %     if debias, ylabel('Normalized debiased MSE', 'fontsize',14); end;
    %     title('EM-BG-GAMP algorithm')
    %     grid on;
    %     legend('1-bit', '2-bit', '3-bit', '4-bit', '\infty-bit')
    
    %     figure,
    %     subplot(121)
    %     gtrue2 = reshape(permute(reshape(gtrue,Nr,Nt,Nd),[1,3,2]),Nd*Nr,Nt);
    %     bar3(abs(gtrue2));
    %     title('True value of the channel', 'Fontsize', 16)
    %     xlabel('Tx angle', 'Fontsize', 14)
    %     ylabel('Rx angle & delay', 'Fontsize', 14)
    %
    %     subplot(122)
    %     gtrue2 = reshape(permute(reshape(ghat_BG,Nr,Nt,Nd),[1,3,2]),Nd*Nr,Nt);
    %     bar3(abs(gtrue2));
    %     title('Estimate of the channel', 'Fontsize', 16)
    %     xlabel('Tx angle', 'Fontsize', 14)
    %     ylabel('Rx angle & delay', 'Fontsize', 14)
    
    
    %     figure,
    %     gtrue2 = reshape(permute(reshape(gtrue,Nr,Nt,Nd),[1,3,2]),Nd*Nr,Nt);
    %     bar3(abs(gtrue2));
    %     xlabel('Tx angle', 'Fontsize', 14)
    %     ylabel('Rx angle & delay', 'Fontsize', 14)
    %
    %     figure,
    %     gtrue2 = reshape(permute(reshape(ghat_GAMP,Nr,Nt,Nd),[1,3,2]),Nd*Nr,Nt);
    %     bar3(abs(gtrue2));
    %     xlabel('Tx angle', 'Fontsize', 14)
    %     ylabel('Rx angle & delay', 'Fontsize', 14)
end;

%% EM-GM-GAMP
if EM_GM_GAMP_flag == 1
    %%%% GAMP
    NMSE_GM_GAMP_dB = 10*log10(mean(MSE_GM_GAMP, 4)/mean(Gnorm.^2));
    % MSE_GAMP_dB = 10*log10(MSE_GAMP(:,:,:,4)./Gnorm(:,:,:,end));
    %     GI_MSE_GM_dB = 10*log10(mean(GI_MSE_GM, 4));
    
    %     % MSE versus SNR
    %     figure,
    %     plot(SNRdB_vector, NMSE_GM_GAMP_dB(:,:, end), '.-');
    %     xlabel('SNR', 'fontsize',14)
    %     ylabel('Normalized MSE', 'fontsize',14)
    %     if debias, ylabel('Normalized debiased MSE', 'fontsize',14); end;
    %     title('EM-GM-GAMP algorithm')
    %     grid on;
    %     legend('1-bit', '2-bit', '3-bit', '4-bit', '\infty-bit')
    
    %     figure,
    %     subplot(121)
    %     gtrue2 = reshape(permute(reshape(gtrue,Nr,Nt,Nd),[1,3,2]),Nd*Nr,Nt);
    %     bar3(abs(gtrue2));
    %     title('True value of the channel', 'Fontsize', 16)
    %     xlabel('Tx angle', 'Fontsize', 14)
    %     ylabel('Rx angle & delay', 'Fontsize', 14)
    %
    %     subplot(122)
    %     gtrue2 = reshape(permute(reshape(ghat_GM,Nr,Nt,Nd),[1,3,2]),Nd*Nr,Nt);
    %     bar3(abs(gtrue2));
    %     title('Estimate of the channel', 'Fontsize', 16)
    %     xlabel('Tx angle', 'Fontsize', 14)
    %     ylabel('Rx angle & delay', 'Fontsize', 14)
    
    
    %     figure,
    %     gtrue2 = reshape(permute(reshape(gtrue,Nr,Nt,Nd),[1,3,2]),Nd*Nr,Nt);
    %     bar3(abs(gtrue2));
    %     xlabel('Tx angle', 'Fontsize', 14)
    %     ax = gca;
    %     ax.XTick = [1 8 16];
    %     ylim([0 Nd*Nr])
    %     ylabel('Rx angle & delay', 'Fontsize', 14)
    %     zlabel('$\left|\mathbf{X}_{i,j} \right|$','interpreter','latex')
    %
    %     figure,
    %     ghat2 = reshape(permute(reshape(ghat_GM,Nr,Nt,Nd),[1,3,2]),Nd*Nr,Nt);
    %     bar3(abs(ghat2));
    %     xlabel('Tx angle', 'Fontsize', 14)
    %     ax = gca;
    %     ax.XTick = [1 8 16];
    %     ylim([0 Nd*Nr])
    %     ylabel('Rx angle & delay', 'Fontsize', 14)
    %     zlabel('$\left|\widehat{\mathbf{X}}_{i,j} \right|$','interpreter','latex')
    %
    %     figure,
    %     bar3(abs(gtrue2 - ghat2));
    %     xlabel('Tx angle', 'Fontsize', 14)
    %     ax = gca;
    %     ax.XTick = [1 8 16];
    %     ylim([0 Nd*Nr])
    %     ylabel('Rx angle & delay', 'Fontsize', 14)
    %     zlabel('$\left|\mathbf{X}_{i,j}-\widehat{\mathbf{X}}_{i,j} \right|$','interpreter','latex')
    %
    %     figure,
    %     bar3(abs(gtrue2 - ghat2)./abs(gtrue2));
    %     xlabel('Tx angle', 'Fontsize', 14)
    %     ax = gca;
    %     ax.XTick = [1 8 16];
    %     ylim([0 Nd*Nr])
    %     ylabel('Rx angle & delay', 'Fontsize', 14)
    %     zlabel('$\frac{\left|\mathbf{X}_{i,j}-\widehat{\mathbf{X}}_{i,j} \right|}{\left|\mathbf{X}_{i,j}\right|}$','interpreter','latex')
    %
    %     figure,
    %     stem(abs(gtrue));
    %     hold on;
    %     stem(abs(ghat_GM));
end;

%% LS
if LS_flag == 1
    NMSE_LS_dB = 10*log10(mean(MSE_LS, 4)/mean(Gnorm.^2));
    % MSE_LS_dB = 10*log10(MSE_LS(:,:,:,4)./Gnorm(:,:,:,end));
    %     MSE_LS_dB = 10*log10(mean(GI_MSE_LS, 4));
    
    %     % MSE versus SNR
    %     figure,
    %     plot(SNRdB_vector, NMSE_LS_dB(:,:, end), '.-');
    %     xlabel('SNR', 'fontsize',14)
    %     ylabel('Normalized MSE', 'fontsize',14)
    %     if debias, ylabel('Normalized debiased MSE', 'fontsize',14); end;
    %     title('LS algorithm')
    %     grid on;
    %     legend('1-bit', '2-bit', '3-bit', '4-bit', '\infty-bit')
    
    %     % MSE versus Np
    %     figure,
    %     semilogx(Np_vector, ( reshape( NMSE_LS_dB(end,:,:), length(bit_vector), []) )');
    %     xlabel('Training Block Length', 'fontsize',14)
    %     ylabel('Normalized MSE', 'fontsize',14)
    %     if debias, ylabel('Normalized debiased MSE', 'fontsize',14); end;
    %     title('LS algorithm')
    %     grid on;
    %
    %     figure,
    %     subplot(121)
    %     gtrue2 = reshape(permute(reshape(gtrue,Nr,Nt,Nd),[1,3,2]),Nd*Nr,Nt);
    %     bar3(abs(gtrue2));
    %     title('True value of the channel', 'Fontsize', 16)
    %     xlabel('Tx angle', 'Fontsize', 14)
    %     ylabel('Rx angle & delay', 'Fontsize', 14)
    %
    %     subplot(122)
    %     gtrue2 = reshape(permute(reshape(ghat_LS,Nr,Nt,Nd),[1,3,2]),Nd*Nr,Nt);
    %     bar3(abs(gtrue2));
    %     title('Estimate of the channel', 'Fontsize', 16)
    %     xlabel('Tx angle', 'Fontsize', 14)
    %     ylabel('Rx angle & delay', 'Fontsize', 14)
end;

%% LMMSE
if LMMSE_flag == 1
    NMSE_LMMSE_dB = 10*log10(mean(MSE_LMMSE, 4)/mean(Gnorm.^2));
    % MSE_LS_dB = 10*log10(MSE_LS(:,:,:,4)./Gnorm(:,:,:,end));
    %     MSE_LS_dB = 10*log10(mean(GI_MSE_LS, 4));
    
    %     % MSE versus SNR
    %     figure,
    %     plot(SNRdB_vector, NMSE_LMMSE_dB(:,:, end), '.-');
    %     xlabel('SNR', 'fontsize',14)
    %     ylabel('Normalized MSE', 'fontsize',14)
    %     if debias, ylabel('Normalized debiased MSE', 'fontsize',14); end;
    %     title('LMMSE algorithm')
    %     grid on;
    %     legend('1-bit', '2-bit', '3-bit', '4-bit', '\infty-bit')
    
    %     % MSE versus Np
    %     figure,
    %     semilogx(Np_vector, ( reshape( NMSE_LMMSE_dB(end,:,:), length(bit_vector), []) )');
    %     xlabel('Training Block Length', 'fontsize',14)
    %     ylabel('Normalized MSE', 'fontsize',14)
    %     if debias, ylabel('Normalized debiased MSE', 'fontsize',14); end;
    %     title('LMMSE algorithm')
    %     grid on;
    %
    %     figure,
    %     subplot(121)
    %     gtrue2 = reshape(permute(reshape(gtrue,Nr,Nt,Nd),[1,3,2]),Nd*Nr,Nt);
    %     bar3(abs(gtrue2));
    %     title('True value of the channel', 'Fontsize', 16)
    %     xlabel('Tx angle', 'Fontsize', 14)
    %     ylabel('Rx angle & delay', 'Fontsize', 14)
    %
    %     subplot(122)
    %     gtrue2 = reshape(permute(reshape(ghat_LMMSE,Nr,Nt,Nd),[1,3,2]),Nd*Nr,Nt);
    %     bar3(abs(gtrue2));
    %     title('Estimate of the channel', 'Fontsize', 16)
    %     xlabel('Tx angle', 'Fontsize', 14)
    %     ylabel('Rx angle & delay', 'Fontsize', 14)
    
end;


%% BPDN
if BPDN_flag == 1
    NMSE_BPDN_dB = 10*log10(mean(MSE_BPDN, 4)/mean(Gnorm.^2));
    % MSE_LS_dB = 10*log10(MSE_LS(:,:,:,4)./Gnorm(:,:,:,end));
    % MSE_LS_dB = 10*log10(mean(GI_MSE_LS, 4));
    
    %     % MSE versus SNR
    %     figure,
    %     plot(SNRdB_vector, NMSE_BPDN_dB(:,:, end), '.-');
    %     xlabel('SNR', 'fontsize',14)
    %     ylabel('Normalized MSE', 'fontsize',14)
    %     if debias, ylabel('Normalized debiased MSE', 'fontsize',14); end;
    %     title('BPDN algorithm')
    %     grid on;
    %     legend('1-bit', '2-bit', '3-bit', '4-bit', '\infty-bit')
end;


%% Lasso
if Lasso_flag == 1
    NMSE_Lasso_dB = 10*log10(mean(MSE_Lasso, 4)/mean(Gnorm.^2));
    
    
    %     % MSE versus SNR
    %     figure,
    %     plot(SNRdB_vector, NMSE_Lasso_dB(:,:, end), '.-');
    %     xlabel('SNR', 'fontsize',14)
    %     ylabel('Normalized MSE', 'fontsize',14)
    %     if debias, ylabel('Normalized debiased MSE', 'fontsize',14); end;
    %     title('Lasso algorithm')
    %     grid on;
    %     legend('1-bit', '2-bit', '3-bit', '4-bit', '\infty-bit')
end;

%% QIST
if QIST_flag == 1
    NMSE_QIST_dB = 10*log10(mean(MSE_QIST, 4)/mean(Gnorm.^2));
    %     % MSE versus SNR
    %     figure,
    %     plot(SNRdB_vector, NMSE_QIST_dB(:,:, end), '.-');
    %     xlabel('SNR', 'fontsize',14)
    %     ylabel('Normalized MSE', 'fontsize',14)
    %     if debias, ylabel('Normalized debiased MSE', 'fontsize',14); end;
    %     title('QIST algorithm')
    %     grid on;
    %     legend('1-bit', '2-bit', '3-bit', '4-bit', '\infty-bit')
end;

%% QIHT
if QIHT_flag == 1
    NMSE_QIHT_dB = 10*log10(mean(MSE_QIHT, 4)/mean(Gnorm.^2));
    %     % MSE versus SNR
    %     figure,
    %     plot(SNRdB_vector, NMSE_QIHT_dB(:,:, end), '.-');
    %     xlabel('SNR', 'fontsize',14)
    %     ylabel('Normalized MSE', 'fontsize',14)
    %     if debias, ylabel('Normalized debiased MSE', 'fontsize',14); end;
    %     title('QIHT algorithm')
    %     grid on;
    %     legend('1-bit', '2-bit', '3-bit', '4-bit', '\infty-bit')
end;


%% CVTOSHT
if CVTOSHT_flag == 1
    NMSE_CVTOSHT_dB = 10*log10(mean(MSE_CVTOSHT, 4)/mean(Gnorm.^2));
    % MSE_CVTOSHT_dB = 10*log10(MSE_CVTOSHT(:,:,:,4)./Gnorm(:,:,:,end));
    % MSE_CVTOSHT_dB = 10*log10(mean(GI_MSE_CVTOSHT, 4));
    %
    %     % MSE versus SNR
    %     figure,
    %     plot(SNRdB_vector, NMSE_CVTOSHT_dB(:,:, end), '.-');
    %     xlabel('SNR', 'fontsize',14)
    %     ylabel('Normalized MSE', 'fontsize',14)
    %     if debias, ylabel('Normalized debiased MSE', 'fontsize',14); end;
    %     title('CVTOSHT algorithm')
    %     grid on;
    %     legend('1-bit', '2-bit', '3-bit', '4-bit', '\infty-bit')
end;

%
% plot(Rate_pcsi);
% hold on;
% plot(Rate_BG);
% plot(Rate_GM);
% plot(Rate_LS);
% plot(Rate_LMMSE);
% plot(Rate_BPDN);

% profsave(profile('info'), 'UPA_64_64_16_2clusters_2048_Toeplitz')

if compute_rate_flag == 1
    MI_pcsi_mean = mean(MI_pcsi,4);
    MI_BG_VAMP_mean = mean(MI_BG_VAMP,4);
    MI_GM_VAMP_mean = mean(MI_GM_VAMP,4);
    MI_BG_GAMP_mean = mean(MI_BG_GAMP,4);
    MI_GM_GAMP_mean = mean(MI_GM_GAMP,4);
    MI_LS_mean = mean(MI_LS,4);
    MI_LMMSE_mean = mean(MI_LMMSE,4);
    MI_BPDN_mean = mean(MI_BPDN,4);
    %     MI_QIST_mean = mean(MI_QIST,4);
    MI_QIHT_mean = mean(MI_QIHT,4);
    
    Rate_pcsi_mean = mean(Rate_pcsi,4);
    Rate_BG_VAMP_mean = mean(Rate_BG_VAMP,4);
    Rate_GM_VAMP_mean = mean(Rate_GM_VAMP,4);
    Rate_BG_GAMP_mean = mean(Rate_BG_GAMP,4);
    Rate_GM_GAMP_mean = mean(Rate_GM_GAMP,4);
    Rate_LS_mean = mean(Rate_LS,4);
    Rate_LMMSE_mean = mean(Rate_LMMSE,4);
    Rate_BPDN_mean = mean(Rate_BPDN,4);
    %     Rate_QIST_mean = mean(Rate_QIST,4);
    Rate_QIHT_mean = mean(Rate_QIHT,4);
    
    figure,
    plot(SNRdB_vector, Rate_pcsi_mean(:,:,end));
    grid on;
    hold on;
    plot(SNRdB_vector, Rate_BG_VAMP_mean(:,:,end), '.-');
    plot(SNRdB_vector, Rate_GM_VAMP_mean(:,:,end), 's-');
    plot(SNRdB_vector, Rate_BG_GAMP_mean(:,:,end), '^-');
    plot(SNRdB_vector, Rate_GM_GAMP_mean(:,:,end), '*-');
    plot(SNRdB_vector, Rate_LS_mean(:,:,end), 'o-');
    plot(SNRdB_vector, Rate_LMMSE_mean(:,:,end), 'v-');
    plot(SNRdB_vector, Rate_BPDN_mean(:,:,end), '.-');
    plot(SNRdB_vector, Rate_QIHT_mean(:,:,end), '.-');
    xlabel('SNR', 'fontsize',14)
    ylabel('Achievable Rate (bps/Hz)', 'fontsize',14)
end;

Results.NMSE.NMSE_BG_VAMP_dB = NMSE_BG_VAMP_dB;
Results.NMSE.NMSE_GM_VAMP_dB = NMSE_GM_VAMP_dB;
Results.NMSE.NMSE_BG_GAMP_dB = NMSE_BG_GAMP_dB;
Results.NMSE.NMSE_GM_GAMP_dB = NMSE_GM_GAMP_dB;
Results.NMSE.NMSE_LS_dB = NMSE_LS_dB;
Results.NMSE.NMSE_LMMSE_dB = NMSE_LMMSE_dB;
Results.NMSE.NMSE_BPDN_dB = NMSE_BPDN_dB;
Results.NMSE.NMSE_QIHT_dB = NMSE_QIHT_dB;
Results.MI.MI_pcsi_mean = MI_pcsi_mean;
Results.MI.MI_BG_VAMP_mean = MI_BG_VAMP_mean;
Results.MI.MI_GM_VAMP_mean = MI_GM_VAMP_mean;
Results.MI.MI_BG_GAMP_mean = MI_BG_GAMP_mean;
Results.MI.MI_GM_GAMP_mean = MI_GM_GAMP_mean;
Results.MI.MI_LS_mean = MI_LS_mean;
Results.MI.MI_LMMSE_mean = MI_LMMSE_mean;
Results.MI.MI_BPDN_mean = MI_BPDN_mean;
Results.MI.MI_QIHT_mean = MI_QIHT_mean;
Results.Rate.Rate_BG_VAMP_mean = Rate_BG_VAMP_mean;
Results.Rate.Rate_GM_VAMP_mean = Rate_GM_VAMP_mean;
Results.Rate.Rate_BG_GAMP_mean = Rate_BG_GAMP_mean;
Results.Rate.Rate_GM_GAMP_mean = Rate_GM_GAMP_mean;
Results.Rate.Rate_LS_mean = Rate_LS_mean;
Results.Rate.Rate_LMMSE_mean = Rate_LMMSE_mean;
Results.Rate.Rate_BPDN_mean = Rate_BPDN_mean;
Results.Rate.Rate_QIHT_mean = Rate_QIHT_mean;
end
