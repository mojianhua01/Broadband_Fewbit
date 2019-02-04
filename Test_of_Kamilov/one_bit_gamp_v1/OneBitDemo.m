% ONEBITDEMO simulates recovery of random variables distributed
% according to Gauss-Bernoulli distribution. Thresholds are set to zero or
% adaptively using the centroids.
%
% Corresponding paper: 
% U. S. Kamilov, A. Bourquard, A. Amini, and M. Unser,
% "One-Bit Measurements with Adaptive Thresholds," 
% IEEE Signal Process. Letters, vol. 19, no. 10, pp. 607-610, 
% October 2012.
%
% Ulugbek S. Kamilov, BIG, EPFL, 2012.

clear; close all; clc;

% Distortion metrics
computeMse = @(noise) 10*log10((norm(noise(:))^2)/length(noise));
computeSnr = @(sig, noise) 10*log10((norm(sig(:))^2)/(norm(noise(:))^2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Random stream
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set random stream so we can reproduce the same experiments5489
%s = RandStream('mt19937ar','Seed', 0);
%RandStream.setGlobalStream(s);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimation Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Length of the signal
n = 100;

% Sparsity ratio
rho = 0.5;

% Signal variance
sigvar = 1;

% Noise variance
v = 0;

% Numer of bits per signal component
betas = linspace(0.1, 5, 5);

% Number of betas
nBetas = length(betas);

% Number of trials
nTrials = 50;

% GAMP iterations
T = 200;

% Initialization for GAMP
xhat0 = 0;
vx0 = rho*sigvar;
init0 = [xhat0; vx0];

% Denoising function
denfunc = @(rhat, vr) denoiseGaussBernoulli(rhat, vr, rho, sigvar);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Store MSEs
msezero = zeros(nTrials, nBetas);
msegamp = zeros(nTrials, nBetas);

% Store SNRs
snrzero = zeros(nTrials, nBetas);
snrgamp = zeros(nTrials, nBetas);

for iTrial = 1:nTrials
    for iBeta = 1:nBetas
        
        % Measurement ratio (beta = m/n) [also nb. of bits per signal element]
        beta = betas(iBeta);
        
        % Number of measurements
        m = round(n*beta);
        
        % Generate the signal (Gaussian with proba. rho)
        x = binornd(1, rho, n, 1) .* randn(n, 1) .* sqrt(sigvar);
        
        % Generate the measurement matrix
        A = (1/sqrt(m)) .* randn(m, n);
        
        % Obtain clean measurements
        z = A*x;
        
        % Noise
        noise = sqrt(v) * randn(m, 1);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Zero-threshold reconstruction
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        xhatzero = signGamp(A, sign(z+noise), zeros(m, 1), init0, denfunc, v, T);
        
        msezero(iTrial, iBeta) = computeMse(xhatzero-x);
        snrzero(iTrial, iBeta) = computeSnr(x, x-xhatzero);
        
        fprintf('[%d/%d][beta: %f bits][Zero-Threshold][MSE = %.4f][SNR = %.4f]\n',...
            iTrial, nTrials, beta, msezero(iTrial, iBeta), snrzero(iTrial, iBeta));        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Adaptive-GAMP thresholds
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Initial reconstruction
        xhatgamp = adaptiveGamp(A, z, init0, denfunc,...
            'verbose', false,...
            'oracle', x,...
            'blockSize', 1,...
            'nIterations', 20,...
            'plotReconstruction', false,...
            'NoiseVariance', v,...,
            'minit', 1,...
            'tauinit', 0,...
            'variableInitialization', false);
        
        msegamp(iTrial, iBeta) = computeMse(xhatgamp-x);
        snrgamp(iTrial, iBeta) = computeSnr(x, x-xhatgamp);
        
        fprintf('[%d/%d][beta: %f bits][Gamp-Threshold][MSE = %.4f][SNR = %.4f]\n',...
            iTrial, nTrials, beta, msegamp(iTrial, iBeta), snrgamp(iTrial, iBeta));
        
        save('OneBitDemo.mat');
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Figure SNR
figure('Color', 'w', 'Name', 'Bits-vs-SNR');

% Plot reconstruction
plot(...
    betas, mean(snrzero, 1), 'sr--',...
    betas, mean(snrgamp, 1), 'ob-',...
    'LineWidth', 1.5);
xlim([min(betas), max(betas)]);
set(gca, 'FontSize', 14);
title(sprintf('Reconstruction SNR: n = %d; rho = %.2f', n, rho));
legend('Zero-Threshold', 'GAMP-Thresholds');
grid on;
ylabel('SNR (dB)');
xlabel('Rate (bits/signal entry)');