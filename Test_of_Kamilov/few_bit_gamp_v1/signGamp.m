function [xhat, vx] = signGamp(A, y, tau, init, denfunc, vn, T, oracle)
% SIGNGAMP reconstruction of a vector from single-bit (sign) measurements
%
% [xhat, vx] = signGamp(A, y, tau, denfunc, vn, T, oracle)
%
% Input:
% - A: measurement matrix (m x n)
% - y: sign measurements (+1 or -1) (m x 1)
% - tau: quantizer thresholds
% - init: initialization for the signal and variance [xhat0; vx0]
% - denfunc: prior dependent denoising function (handle)
% - vn: AWGN variance
% - T: number of iterations
% - oracle: original signal
%
% Output:
% - xhat: reconstructed signal (n x 1)
% - vx: predicted MSE (n x 1)
%
% Ulugbek Kamilov, BIG, EPFL, 2012.

% Number of measurements and dimension of the signal
[m, n] = size(A);

% MSE
computeMse = @(noise) 10*log10((norm(noise(:))^2)/n);
computeSnr = @(sig, noise) 10*log10((norm(sig(:))^2)/(norm(noise(:))^2));
absdiff = @(x) sum(abs(x))/length(x);

% Initialize the estimates
if(numel(init) == 2)
    xhat = init(1)*zeros(n, 1);
    vx = init(2)*ones(n, 1);
else
    xhat = init(:, 1);
    vx = init(:, 2);
end

% Initialize shat
shat = zeros(m, 1);

dampFac = 0.5;

% Previous estimate
xhatprev = xhat;
shatprev = shat;
vxprev = vx;

% Number of unchanged iterations
count = 3;

% Hadamard product of the matrix
AA = A.*A;

% If we have access to the oracle
if(exist('oracle', 'var'))
    % Create new figure
    h = figure;
end

% Perform estimation
for t = 1:T
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Measurement update
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Linear
    vp = AA*vx;
    phat = A*xhat - vp.*shat;
    
    % Truncated Gaussian
    [ez, vz] = truncatedGaussianMoments(y, tau+phat, vp+vn);
    
    % Non-Linear
    shat = (ez - (tau+phat))./(vp+vn);
    vs = (1 - vz./(vp+vn))./(vp+vn);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Estimation update
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Linear
    vr = 1./(AA' * vs);
    rhat = xhat + vr .* (A' * shat);
    
    % Non-linear
    [xhat, vx] = denfunc(rhat, vr);
    
    %Damp
    xhat = dampFac*xhat + (1-dampFac)*xhatprev;
    shat = dampFac*shat + (1-dampFac)*shatprev;
    vx = dampFac*vx + (1-dampFac)*vxprev;
    
    % If without a change
    if(absdiff(xhat - xhatprev) < 1e-6)
        count = count - 1;
    end
    
    % Save previous xhat
    xhatprev = xhat;
    shatprev = shat;
    vxprev = vx;
    
    if(exist('oracle', 'var'))
        mse = computeMse(oracle-xhat);
        snr = computeSnr(oracle, oracle-xhat);
        
        % Print messages
        fprintf('[%d/%d][mse = %.4f][snr = %.4f]\n', t, T, mse, snr);
        
        % Plot
        figure(h);
        subplot(3, 1, 1);
        plot(t, mse, '.');
        hold on;
        title(sprintf('MSE: %.4f', mse));
        xlim([1 T]);
        
        subplot(3, 1, 2);
        plot(t, snr, '.');
        hold on;
        title(sprintf('SNR: %.4f', snr));
        xlim([1 T]);
        
        subplot(3, 1, 3);
        plot(1:n, oracle, 'o', 1:n, xhat, '.');
        title('Estimate');
        legend('Oracle', 'xhat');
        xlim([1 n]);
        drawnow;
    end
    
    % Stopping criterion
    if(count <= 0)
        break;
    end
    
end