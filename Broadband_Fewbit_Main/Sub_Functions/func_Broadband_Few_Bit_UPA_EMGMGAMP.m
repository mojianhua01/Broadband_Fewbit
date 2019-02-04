function [ ghat, time_EM_GM_GAMP] = func_Broadband_Few_Bit_UPA_EMGMGAMP( T_cyclic, r_bit, Params)
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

Gtrue = Params.Gtrue;

% Bt = Params.Bt;
% V_cyclic_0 = kron(eye(Nd), Bt') * T_cyclic;
% fast implementation of V_cyclic = kron(eye(Nd), Bt') * T_cyclic
T_temp = reshape(T_cyclic, Nt, Nd, Np);
V_temp = Bhfast(T_temp, Nt);
V_cyclic = reshape(V_temp, Nd*Nt, Np);

% EG-BG-GAMP setup
sparseRat = 0.05;    % Bernoulli-Gaussian active probability
BGmean = 0;        	% Bernoulli-Gaussian active mean
BGvar = 1/sparseRat/Nd;         % Bernoulli-Gaussian active variance

% Decide on MAP or MMSE GAMP. Only support MMSE GAMP in our codes!
map = false;

% Auto tune the mean, variance and sparseRat.
% This is extremely important for the algorithm to work!!!

% EG-GM-GAMP setup
if 1
    Lgm = 3;            % number of GM components
    omega0 = ones(Nt*Nr*Nd, 1, Lgm)/Lgm;
    %theta0 = bsxfun(@times, ones(Nt*Nr*Nd, 1, Lgm),...
    %                reshape(linspace(-1,1,Lgm)/sqrt(1/3),[1,1,Lgm]) );
    %phi0 = ones(Nt*Nr*Nd, 1, Lgm)/Lgm;
    theta0 = BGmean*ones(Nt*Nr*Nd, 1, Lgm);
    vars = [1:Lgm]/sqrt(Lgm);
    phi0 = bsxfun(@times, ones(Nt*Nr*Nd, 1, Lgm),...
        reshape(BGvar*vars/mean(vars),[1,1,Lgm]) );
    inputEst0 = CGMEstimIn(omega0, theta0, phi0, 'autoTune', true, 'thetaTune', false);
else % BG case for debugging
    inputEst0 = CAwgnEstimIn(BGmean, BGvar, map, 'autoTune', true, 'mean0Tune', false);
end
inputEst = SparseScaEstim(inputEst0, sparseRat, 0, 'autoTune', true);


if bit == +inf
    outputEst = CAwgnEstimOut(r_bit, wvar, map);
else % uniform quantization
    % Create an output estimation class corresponding to the Gaussian noise.
    % Note that the observation vector is passed to the class constructor, not
    % the gampEst function.
    % outputEst = CQuantizeEstimOut(r_bit, bit, stepsize, 0, wvar, map);
    outputEst = CQuantizeEstimOut(r_bit, bit, stepsize, dither_mean, wvar, map);
end;

% Set the default options
opt = GampOpt();
if bit == +inf
    opt.nit = 50;
    opt.tol = max(min(1e-3, 1/SNR_linear),1e-15);
    opt.stepMax = 0.5;
    opt.step = 0.5;
else
    opt.nit = 25; %25
    opt.tol = 10^(-2.5); % 10^(-2.5)
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

% % stop function: stop when the NMSE is smaller than that of BPDN (i.e., SPGL1) algorithm
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
    Reps = 1;  % 20
    for ii=1:1:Reps
        [estFin] = gampEstBasic(inputEst, outputEst, lA, opt);
    end;
    time_EM_GM_GAMP = toc/Reps
else % debug mode
    tic;
    Reps = 1;
    for i=1:1:Reps
        [estFin] = gampEstBasic(inputEst, outputEst, lA, opt);
    end;
    time_EM_BG_GAMP = toc/Reps
    %     gampShowHist(estHist,optFin,Gtrue(:));
    %     drawnow
    %     clf
    figure,
    plot(10*log10(sum(abs(bsxfun(@minus, estHist.xhat, Gtrue(:))).^2)./norm(Gtrue(:))^2))
    10*log10(sum(abs(bsxfun(@minus, estHist.xhat, Gtrue(:))).^2)/norm(Gtrue(:))^2);
    xlabel('Iteration')
    ylabel('Normalized MSE')
    title('EM-GM-GAMP')
    grid on;
end


ghat = estFin.xhat(:, end);

% Ghat = reshape(ghat, Nr, Nt, Nd);
end

