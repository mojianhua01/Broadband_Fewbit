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

% Extract
minp = muin(ip);
sinp = sqrt(vain(ip));

% Mean to stddev ratio
alpha = -minp ./ sinp;

%lambda = normpdf(alpha) ./ (1-normcdf(alpha));
lambda = normpdf(alpha) ./ (0.5*erfc(alpha/sqrt(2)));
delta = lambda .* (lambda - alpha);

% Moments
e1 = minp + sinp .* lambda;
v1 = (sinp.^2).*(1 - delta);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Negative moments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Extract
minn = muin(in);
sinn = sqrt(vain(in));

% Mean to stddev ratio
beta = -minn ./ sinn;

lambda = normpdf(beta) ./ normcdf(beta);

% Moments
e2 = minn - sinn .* lambda;
v2 = (sinn.^2).*(1 - beta.*lambda - (lambda.^2));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Combine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mout(ip) = e1;
mout(in) = e2;
vout(ip) = v1;
vout(in) = v2;