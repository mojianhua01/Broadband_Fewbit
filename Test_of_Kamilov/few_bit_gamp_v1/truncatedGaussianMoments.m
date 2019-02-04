function [mout, vout] = truncatedGaussianMoments(y, muin, vain)
% TRUNCATEDGAUSSIANMOMENTS returns posterior mean and variance for
% estimation: Y ~ Q(Z) Z~N(muin, vain), with Q(.) 1-bit quantizer.
%
% Input:
% - y: sign measurements
% - muin: prior mean
% - vain:  prior variance
%
% Output:
% - mout: E(Z | Y = y)
% - vout: Var(Z | Y = y)
%
% Ulugbek Kamilov, BIG, EPFL, 2012.

% Split to positive and negative parts
ip = logical(y >= 0);
in = ~ip;

% Initialize
mout = zeros(length(muin), 1);
vout = zeros(length(vain), 1);

% Check for error
%assert(max(abs(muin./sqrt(vain))) <= 7, 'TRUNCATEDGAUSSIANMOMENTS: max(muin./sqrt(vain)) <= 7')
%il = abs(muin./sqrt(vain)) > 7;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Positive moments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
wvar = 1;
% Extract
phat = muin(ip);
pvar = vain(ip);

% Mean to stddev ratio
alpha = -phat ./ sqrt(pvar + wvar);

%lambda = normpdf(alpha) ./ (1-normcdf(alpha));
lambda = normpdf(alpha) ./ (0.5*erfc(alpha/sqrt(2)));

% Moments
e1 = phat + pvar./sqrt(wvar + pvar) .* lambda;
v1 = pvar - (pvar.^2 ./ (pvar + wvar)) .* lambda .*(- alpha + lambda);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Negative moments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Extract
phat = muin(in);
pvar = sqrt(vain(in));

% Mean to stddev ratio
beta = -phat ./ sqrt(pvar + wvar);

lambda = normpdf(beta) ./ normcdf(beta);

% Moments
e2 = phat + pvar./sqrt(wvar + pvar) .* lambda;
v2 = pvar - (pvar.^2 ./ (pvar + wvar)) .* lambda .*( beta + lambda);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Combine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mout(ip) = e1;
mout(in) = e2;
vout(ip) = v1;
vout(in) = v2;