function [ Ghat ] = func_Broadband_Few_Bit_Hybrid( Gtrue, SNR, S_train_hybrid, bit, Analog_combining_pattern)
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
N = size(S_train_hybrid, 2);

% EG-BG-GAMP setup
sparseRat = 0.3;    % Bernoulli-Gaussian active probability
BGmean = 0;        	% Bernoulli-Gaussian active mean
BGvar = 1;         	% Bernoulli-Gaussian active variance

SNR_linear = 10^(SNR/10);

S_train_scaled = S_train_hybrid./norm(S_train_hybrid, 'fro')*sqrt(N)*sqrt(SNR_linear);
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

z_digital = dftmtx(Nr)/sqrt(Nr) * reshape(Gtrue, Nr, []) * S;

z = reshape( z_digital .* Analog_combining_pattern, [], 1);

wvar = 1;

% Compute training class labels based on synthetic "true" hyperplane
y = z + sqrt(1/2)*wvar*randn(Nr*N,2)*[1;1i];

% Decide on MAP or MMSE GAMP. Only support MMSE GAMP in our codes!
map = false;

% Auto tune the mean, variance and sparseRat.
% This is extremely important for the algorithm to work!!!
inputEst = CAwgnEstimIn(BGmean, BGvar, map, 'autoTune', true);
inputEst = SparseScaEstim(inputEst, sparseRat, 0, 'autoTune', true);

stepsize_vector = [1.5956 0.9957 0.586 0.3352 0.1881 0.1041 0.0569 0.0308];

if bit == +inf
    r_bit = y;
    outputEst = CAwgnEstimOut(r_bit, wvar, map);
else
    if bit == 1;
        stepsize = 1.01 * max(abs([real(y); imag(y)]));  % if ADC has only one bit, set the stepsize correctly
    else
        stepsize = ( 1.01 * max(abs([real(y); imag(y)])))/(2^(bit-1));
    end;
    
    %     if bit == 1;
    %         stepsize = max(abs([real(y); imag(y)])) * 1.01;  % if ADC has only one bit, set the stepsize correctly
    %     else
    %         stepsize = rms([real(y); imag(y)])* stepsize_vector(bit);
    %     end;
    
    r_bit_real = min ( max( floor( real(y)/stepsize) + 2^(bit - 1), 0), 2^bit-1) ;
    r_bit_imag = min( max( floor( imag(y)/stepsize )+ 2^(bit - 1), 0), 2^bit-1);
    
    r_bit = r_bit_real + 1j * r_bit_imag;
    
    % Create an output estimation class corresponding to the Gaussian noise.
    % Note that the observation vector is passed to the class constructor, not
    % the gampEst function.
    outputEst = CQuantizeEstimOut(r_bit, bit, stepsize, 0, wvar, map);
end;

% Set the default options
opt = GampOpt();
opt.nit = 40;
opt.tol = 1e-4;
opt.uniformVariance=0;
opt.pvarMin = 1e-5;     % At very high SNR, use very small pvarMin!
opt.xvarMin = 1e-5;     % At very high SNR, use very small pvarMin!
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

hA = @(in) Afast_Hybrid(in,Analog_combining_pattern, S, Nr,1);
hAt = @(in) Afast_Hybrid(in,Analog_combining_pattern, S, Nr,0);
lA = FxnhandleLinTrans(Nr*N, Nt*Nr*Nd, hA, hAt); % linear transform object

% Run the GAMP algorithm
tic
[estFin, ~, estHist] = gampEst(inputEst, outputEst, lA, opt);
timeGAMP = toc

ghat = estFin.xhat(:, end);

%%Plot a figure showing the convergence speed
% for ii=1:1:size(estHist.xhat, 2)
%     Error(ii) = ( norm( estHist.xhat(:, ii) - Gtrue(:)) )^2 / norm(Gtrue(:))^2;
% end;
% figure, 
% plot(Error, 'o-')
% xlabel('Iteration number')
% ylabel('Estimation Error: $\frac{|\widehat{X} - X|^2}{|X|^2}$','Interpreter','LaTex')

Ghat = reshape(ghat, Nr, Nt, Nd);
end

