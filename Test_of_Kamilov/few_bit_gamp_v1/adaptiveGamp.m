function [xhat, vx, tau] = adaptiveGamp(A, z, init, denfunc, varargin)
% ADAPTIVEGAMP simulates a measurement and reconstruction system based on
% GAMP and single-bit quantization with variable thresholds.
%
% [xhat, vx, tau] = adaptiveGamp(A, z, init, denfunc, varargin)
%
% Input:
% - A: measurement matrix (m x n)
% - z: clean measurements z = A*x (m x 1)
% - init: initialization for the signal and variance [xhat0; vx0]
% - denfunc: prior dependent denoising function (handle)
% - varargin: options for the algorithm
%       Options:
%       - Oracle: original signal x (n x 1)
%       - nIterations: number of iterations for GAMP (def: 20)
%       - NoiseVariance: AWGN variance (def: 0)
%       - Verbose: intermediate command-line messages (def: false)
%       - VariableInitializaion: initialize GAMP with previous xhat (def:
%         false)
%       - PlotReconstruction: Display each GAMP reconstruction process
%         (def: false)
%       - Minit: size of initial reconstruction with tau = 0 (def: n)
%       - BlockSize: Number of adaptive measurements per acquisition (def:
%         1)
%       - TauInit: thresholds for initial reconstruction (def: zeros(m, 1))
%
% Output:
% - xhat: final reconstruction (n x 1)
% - vx: posterior variance (n x 1)
% - tau: thresholds (m x 1)
%
% Ulugbek Kamilov, BIG, EPFL, 2012.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Useful shortcuts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

computeMse = @(noise) 10*log10((norm(noise(:))^2)/length(noise));
computeSnr = @(sig, noise) 10*log10((norm(sig(:))^2)/(norm(noise(:))^2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set default values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Problem dimensions
[m, n] = size(A);

% Number of GAMP iterations
nIterations = 20;

% AWGN variance
noiseVariance = 0;

% Initial reconstruction
if(m >= n)
    minit = n;
else
    minit = m;
end

% Initial thresholds
tauinit = zeros(minit, 1);

% Measurements block size
blockSize = 1;

% Re-initialize with reconstruction with each block
variableInitialization = false;

% Plot intermediate reconstructions
plotReconstruction = false;

% Display intermediate messages
verbose = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parse the input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of options
nargs = length(varargin);

% Go through options
for i = 1:2:nargs
    % Extract name/value pair
    name = lower(varargin{i});
    value = varargin{i+1};
    
    switch(name)
        case 'niterations'
            nIterations = value;
        case 'oracle'
            oracle = value;
        case 'noisevariance'
            noiseVariance = value;
        case 'minit'
            minit = value;
        case 'tauinit'
            tauinit = value;
        case 'blocksize'
            blockSize = value;
        case 'variableinitialization'
            variableInitialization = value;
        case 'verbose'
            verbose = value;
        case 'plotreconstruction'
            plotReconstruction = value;
        otherwise
            error('adaptiveGamp: input is not recognized!');
    end
end

% Control parameters
assert(length(tauinit) == minit, 'ADAPTIVEGAMP: TauInit is not set.');
assert(noiseVariance >= 0, 'ADAPTIVEGAMP: noiseVariance >= 0.');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial reconstruction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create noise vector
noise = sqrt(noiseVariance) * randn(m, 1);

% Quantized measurements
y = quantize(z(1:minit)+noise(1:minit), tauinit);

% Initial reconstruction
if(~plotReconstruction)
    [xhat, vx] = signGamp(...
        A(1:minit, :), y, tauinit, init, denfunc, noiseVariance, nIterations);
else
    [xhat, vx] = signGamp(...
        A(1:minit, :), y, tauinit, init, denfunc, noiseVariance, nIterations, oracle);
end

% Print message to command line
if(verbose)
    fprintf('[1:%d / %d]', minit, m);
    if(exist('oracle', 'var'))
        fprintf('[mse = %.4f][snr = %.4f]',...
            computeMse(oracle-xhat), computeSnr(oracle, oracle-xhat));
    end
    fprintf('\n');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adaptively acquire
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Current measurement
mlow = minit + 1;

% Highest measurement in the block
mhigh = mlow + blockSize - 1;

% Initial tau
tau = tauinit;

while(mhigh <= m)
    
    % Thresholds
    tau = [tau; (-A(mlow:mhigh, :)*xhat)];
    
    % Measurements
    y = quantize(z(1:mhigh) + noise(1:mhigh), tau);
    
    assert(min(abs(y))~=0, 'min(abs(y))~=0');
    
    % Initialization
    if(variableInitialization)
        init = [xhat, vx];
    end
    
    % Reconstruct
    if(~plotReconstruction)
        [xhat, vx] = signGamp(...
            A(1:mhigh, :), y, tau, init, denfunc, noiseVariance, nIterations);
    else
        [xhat, vx] = signGamp(...
            A(1:mhigh, :), y, tau, init, denfunc, noiseVariance, nIterations, oracle);
    end
    
    % Print message to command line
    if(verbose)
        fprintf('[%d:%d / %d]', mlow, mhigh, m);
        if(exist('oracle', 'var'))
            fprintf('[mse = %.4f][snr = %.4f]',...
                computeMse(oracle-xhat), computeSnr(oracle, oracle-xhat));
        end
        fprintf('\n');
    end
    
    % Update next block indices
    mlow = mhigh + 1;
    mhigh = mlow + blockSize - 1;
    
end

if(mlow <= m)
    % Thresholds
    tau = [tau; (-A(mlow:m, :)*xhat)];
    
    % Measurements
    y = quantize(z + noise, tau);
    
    assert(min(abs(y))~=0, 'min(abs(y))~=0');
    
    % Initialization
    if(variableInitialization)
        init = [xhat, vx];
    end
    
    % Reconstruct
    if(~plotReconstruction)
        [xhat, vx] = signGamp(...
            A, y, tau, init, denfunc, noiseVariance, nIterations);
    else
        [xhat, vx] = signGamp(...
            A, y, tau, init, denfunc, noiseVariance, nIterations, oracle);
    end
    
    % Print message to command line
    if(verbose)
        fprintf('[%d:%d / %d]', mlow, m, m);
        if(exist('oracle', 'var'))
            fprintf('[mse = %.4f][snr = %.4f]',...
                computeMse(oracle-xhat), computeSnr(oracle, oracle-xhat));
        end
        fprintf('\n');
    end
    
end
