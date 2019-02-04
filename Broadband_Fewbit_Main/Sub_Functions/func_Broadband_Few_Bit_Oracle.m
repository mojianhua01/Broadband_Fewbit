function [ Ghat ] = func_Broadband_Few_Bit( Gtrue, SNR, S_train, bit, Nra, Nre)
% Estimate the channel given the SNR and training signals
% Inputs:
%   S_train:    Nr by N training signals
%   Gtrue:      true value of the channel in the angular domain
%   Bit:        Quantization resolution
%   Stepsize:   Quantization stepsize
%   SNR:        Transmission power in dB
%
% Output:
%   Ghat:       estimate of the channel in the angular domain

% Channel model
Nt = size(Gtrue, 2);
Nr = size(Gtrue, 1);
Nd = size(Gtrue, 3);
N = size(S_train, 2);

% EG-BG-GAMP setup
sparseRat = 0.1;    % Bernoulli-Gaussian active probability
BGmean = 0;        	% Bernoulli-Gaussian active mean
BGvar = 1/sparseRat/Nd;         % Bernoulli-Gaussian active variance

SNR_linear = 10^(SNR/10);

S_train_scaled = S_train./(norm(S_train, 'fro'))*sqrt(N)*sqrt(SNR_linear);
% scale S_train such that sum((S_train_scaled(:,ii)).^2) = SNR_linear for
% ii=1,2,...

% create circular-delay matrices: J{d} rotates down by (d-1) places
J = cell(1,Nd-1);
J{1} = speye(N);
for d=1:Nd-1
    J{d+1} = [spalloc(N-d,d,0), speye(N-d); speye(d), spalloc(d,N-d,0)];
end

S = zeros(Nt*Nd, N);
S(1: Nt,:) = S_train_scaled;
for d=2:1:Nd
    S( (d-1)*Nt+1: d*Nt, :) = S_train_scaled * J{d};
end;

% A = kron(transpose(S), dftmtx(Nr)/sqrt(Nr));
% z = A *gtrue;      % "True" transform vector

z = reshape( dftmtx(Nr)/sqrt(Nr) * reshape(Gtrue, Nr, []) * S, [], 1);

wvar = 1;

% Compute training class labels based on synthetic "true" hyperplane
y = z + sqrt(1/2)*wvar*randn(Nr*N,2)*[1;1i];

% Decide on MAP or MMSE GAMP. Only support MMSE GAMP in our codes!
map = false;


% Auto tune the mean, variance and sparseRat.
% This is extremely important for the algorithm to work!!!

inputEst = CAwgnEstimIn(BGmean, BGvar, map, 'autoTune', false, 'mean0Tune', true);
inputEst = SparseScaEstim(inputEst, sparseRat, 0, 'autoTune', true);
% inputEst = CGMEstimIn(1, 0, 1);
% inputEst = CBGZeroMeanEstimIn(BGvar, sparseRat);

stepsize_vector = [1.5956 0.9957 0.586 0.3352 0.1881 0.1041 0.0569 0.0308];

if bit == +inf
    r_bit = y;
    outputEst = CAwgnEstimOut(r_bit, wvar, map);
else
    if bit == 1;
        stepsize = 1.01 * max(abs([real(y); imag(y)]));  % if ADC has only one bit, set the stepsize correctly
        dither_mean = 0;
        %dither_mean =  rms([real(y); imag(y)]);
        r_bit_real = ( sign(real(y) - dither_mean) + 1 )/2; % positive 1; negative 0
        r_bit_imag = ( sign(imag(y) - dither_mean) + 1 )/2; % 
        r_bit = r_bit_real + 1j * r_bit_imag;
        % Create an output estimation class corresponding to the Gaussian noise.
        % Note that the observation vector is passed to the class constructor, not
        % the gampEst function.
        outputEst = CQuantizeEstimOut(r_bit, bit, stepsize, dither_mean, wvar, map);
%         outputEst = CProbitEstimOut(r_bit, dither_mean, wvar, map);
        
    else
%         stepsize = ( 1.01 * max(abs([real(y); imag(y)])))/(2^(bit-1));
        stepsize = rms([real(y); imag(y)])* stepsize_vector(bit);
        dither_mean = 0;
        
        r_bit_real = min ( max( floor( (real(y) - dither_mean) /stepsize) + 2^(bit - 1), 0), 2^bit-1) ;
        r_bit_imag = min ( max( floor( (imag(y) - dither_mean) /stepsize )+ 2^(bit - 1), 0), 2^bit-1);
        r_bit = r_bit_real + 1j * r_bit_imag;
        
        % Create an output estimation class corresponding to the Gaussian noise.
        % Note that the observation vector is passed to the class constructor, not
        % the gampEst function.
        % outputEst = CQuantizeEstimOut(r_bit, bit, stepsize, 0, wvar, map);
        outputEst = CQuantizeEstimOut(r_bit, bit, stepsize, dither_mean, wvar, map);
    end;

end;

% Set the default options
opt = GampOpt();
opt.nit = 15;
opt.tol = 1e-4; %max(min(1e-3, 10^(-SNR/10)),1e-15);
opt.uniformVariance=0;
% opt.pvarMin = 1e-5;     % At very high SNR, use very small pvarMin!
% opt.xvarMin = 1e-5;     % At very high SNR, use very small pvarMin!
opt.adaptStep = true;
opt.adaptStepBethe= false;
opt.legacyOut=false;
opt.verbose=false;
opt.stepMax=0.25;
opt.step=0.25;

% Demonstrate automatic selection of xvar0
% if 1
%   opt.xvar0auto = true;
%   opt.xhat0 = x + 1*(randn(nx,2)*[1;1j])*norm(x)/sqrt(nx);
% end;

hA = @(in) Afast_Digital_UPA(in,S,Nra, Nre, 1);
hAt = @(in) Afast_Digital_UPA(in,S,Nra, Nre, 0);
lA = FxnhandleLinTrans(Nr*N, Nt*Nr*Nd, hA, hAt); % linear transform object

% Run the GAMP algorithm
tic;
[estFin] = gampEst(inputEst, outputEst, lA, opt);
timeGAMP = toc


ghat = estFin.xhat(:, end);

Ghat = reshape(ghat, Nr, Nt, Nd);
end

_