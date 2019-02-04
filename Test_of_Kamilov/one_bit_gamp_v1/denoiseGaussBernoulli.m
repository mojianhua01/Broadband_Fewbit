function [mout, vout] = denoiseGaussBernoulli(rhat, vr, rho, vx)
% DENOISEGAUSSBERNOULLI returns posterior mean and variance for estimation
% of a GaussBernoulli vector in AWGN.
%
% [mout, vout] = denoiseGaussBernoulli(rhat, vr, rho, vx)
%
% Input:
% - rhat: noisy signal (Rhat = X + N)
% - vr: AWGN variance (Var(N) = vr)
% - rho: probability of zero of the signal (Prob(X = 0) = rho)
% - vx: variance of the Gaussian variable (Var(X | X ~= 0))
%
% Output:
% - mout: E(X | Rhat = rhat)
% - vout: Var(X | Rhat = rhat)
%        where Rhat = X + N, with N ~ N(0, vr)
%
% Ulugbek Kamilov, BIG, EPFL, 2012.

% Variance ratio
rat = vx ./ vr;

% Gaussian ratio
gr = exp(-0.5 .* rat .* (rhat.^2) ./ (vr+vx) + 0.5*log(1 + rat));

% Temporary
temp = (rho./(rho + (1-rho) .* gr)) .* (vx ./ (vx+vr));

% Compute E{X|Q=q}
mout = temp .* rhat;

% Compute Var{X|Q=q}
vout = temp .* (vr + (vx ./ (vx+vr)).*(rhat.^2) - temp.*(rhat.^2));