function [ ghat,  time_EM_BG_GAMP] = func_Broadband_Few_Bit_UPA_EMBGGAMP( T_cyclic, r_bit, Params)
% Estimate the channel given the training signals and quantization outputs
% Inputs:
%   T_cyclic:   training signals
%   r_bit:      Quantization resolution
%   Params:     Parameters
% Output:
%   ghat:       estimate of the channel in the angular domain

% Channel model
Nt = Params.Nt;
Nr = Params.Nr;
Nd = Params.Nd;
Np = Params.Np;

bit = Params.bit;
stepsize = Params.stepsize;
SNR_linear = Params.SNR_linear;
dither_mean = Params.dither_mean;
wvar = Params.wvar;

Bt = Params.Bt;
Gtrue = Params.Gtrue;

% V_cyclic_0 = kron(eye(Nd), Bt') * T_cyclic;

% % 1st fast implementation of V_cyclic = kron(eye(Nd), Bt') * T_cyclic
%                 T_temp = reshape(T_cyclic, Nt, Nd, Np);
%                 V_temp = Bhfast(T_temp, Nt);
%                 V_cyclic_1 = reshape(V_temp, Nd*Nt, Np);

% % 2nd fast implementation of V_cyclic = kron(eye(Nd),
% Bt') * T_cyclic. The time consumption of 1st and 2nd
% implementation are similar.
% V_cyclic = NaN(Nt*Nd, Np);
% V_train = Bhfast(T_cyclic(1:Nt,:), Nt);
V_train = Bt'*T_cyclic(1:Nt,:);
for ii=0:1:Nd-1
    V_cyclic( ii*Nt+1: (ii+1)*Nt, :) = circshift(V_train,ii,2);
end

% we cannot re-order V_cyclic's rows to make it become a Toeplitz matrix

% EG-BG-GAMP setup
sparseRat = 0.05;    % Bernoulli-Gaussian active probability
BGmean = 0;        	% Bernoulli-Gaussian active mean
BGvar = 1/sparseRat/Nd;         % Bernoulli-Gaussian active variance

% Decide on MAP or MMSE GAMP. Only support MMSE GAMP in our codes!
map = false;

% Auto tune the mean, variance and sparseRat.
% This is extremely important for the algorithm to work!!!

% inputEst0 = CAwgnEstimIn(BGmean, BGvar, map, 'autoTune', false, 'mean0Tune', true);
inputEst0 = CAwgnEstimIn(BGmean, BGvar, map, 'autoTune', true, 'mean0Tune', false);
inputEst = SparseScaEstim(inputEst0, sparseRat, 0, 'autoTune', true);
% inputEst = CGMEstimIn(1, 0, 1);
% inputEst = CBGZeroMeanEstimIn(BGvar, sparseRat);

if bit == +inf
    outputEst = CAwgnEstimOut(r_bit, wvar, map);
else
    % Create an output estimation class corresponding to the Gaussian noise.
    % Note that the observation vector is passed to the class constructor, not
    % the gampEst function.
    outputEst = CQuantizeEstimOut(r_bit, bit, stepsize, dither_mean, wvar, map);
    %         outputEst = CProbitEstimOut(r_bit, dither_mean, wvar, map);
end;

% Set the default options
opt = GampOpt();
if bit == +inf
    opt.nit = 50;
    opt.tol = max(min(1e-3, 1/SNR_linear),1e-15);
    opt.stepMax = 0.5;
    opt.step = 0.5;
else
    opt.nit = 25; %50
    opt.tol = 10^(-2.5);  %10^(-10)
    opt.stepMax = 0.5;
    opt.step = 0.5;
end
opt.uniformVariance=1;
% opt.pvarMin = 1e-5;     % At very high SNR, use very small pvarMin!
% opt.xvarMin = 1e-5;     % At very high SNR, use very small pvarMin!
opt.adaptStep = false;
opt.adaptStepBethe= false;
opt.legacyOut=false;
opt.verbose=false;

%% stop function: stop when the NMSE is smaller than that of BPDN (i.e., SPGL1) algorithm
% NMSE_BPDN_dB_vector=[-11.67099606	-12.57338142	-12.91700443	-13.77483569	...
%     -14.29506195	-14.66714104	-14.95196809	-15.34087005	-15.5699147]; % 0dB 4-bit shifted-ZC
% Np_vector = 512*(2:1:10);
% NMSE_BPDN_dB = NMSE_BPDN_dB_vector( Np_vector==Np );
% opt.stopFcn = @(val, xhatFinal, xhatPrevFinal, AxhatFinal)...
%     (10*log10( sum(abs(xhatFinal-Gtrue(:)).^2,1)./sum(abs(Gtrue(:)).^2,1) )<NMSE_BPDN_dB);


%% fast algorithm

%% UPA
hA = @(in) Afast_Digital_UPA(in, T_cyclic, V_cyclic, Params, 1);
hAt = @(in) Afast_Digital_UPA(in, T_cyclic, V_cyclic, Params, 0);

%% ULA
% hA = @(in) Afast_Digital_ULA(in,S,Nr, 1);
% hAt = @(in) Afast_Digital_ULA(in,S,Nr, 0);

lA = FxnhandleLinTrans(Nr*Np, Nt*Nr*Nd, hA, hAt); % linear transform object

% Run the GAMP algorithm
if 1
    tic;
    Reps = 1; %20
    for i=1:1:Reps
        [estFin] = gampEstBasic(inputEst, outputEst, lA, opt);
    end;
    time_EM_BG_GAMP = toc/Reps
else % debug mode
    tic;
    Reps = 1;
    for i=1:1:Reps
        [estFin,optFin,estHist] = gampEst(inputEst, outputEst, lA, opt);
    end;
    time_EM_BG_GAMP = toc/Reps
%         gampShowHist(estHist,optFin,Gtrue(:));
%         drawnow
        clf
        figure,
        plot(sum(abs(bsxfun(@minus, estHist.xhat, Gtrue(:))).^2)./norm(Gtrue(:))^2)
        10*log10( sum(abs(bsxfun(@minus, estHist.xhat, Gtrue(:))).^2)/norm(Gtrue(:))^2)
        xlabel('Iteration')
        ylabel('Normalized MSE')
        title('EM-BG-GAMP')
        grid on;
end


ghat = estFin.xhat(:, end);

% Ghat = reshape(ghat, Nr, Nt, Nd);
end

