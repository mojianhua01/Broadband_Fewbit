function y = quantize(z, tau)
% QUANTIZE single-bit measurements with variable thresholds
%
% y = quantize(z)
% y = quantize(z, tau)
%
% Input:
% - z: array to quantize (m x 1)
% - tau: thresholds (m x 1) [Default: zeros(m, 1)]
%
% Output:
% - y = sign(z + tau)
%
% Ulugbek Kamilov, BIG, EPFL, 2012.

% If threshold was not set
if(~exist('tau', 'var'))    
    tau = 0;
end

% Quantize
y = sign(z + tau);