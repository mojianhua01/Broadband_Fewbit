function [ ghat_ALMMSE ] = func_Broadband_Few_Bit_UPA_ALMMSE( T_cyclic, r_bit, Params)
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
beta = Params.beta;
SNR_linear = Params.SNR_linear;
dither_mean = Params.dither_mean;
wvar = Params.wvar;

% beta_vector = [0.3634 0.1188 0.03744 0.01154];
% if bit == +inf
%     beta = 0;
% else
%     beta = beta_vector(bit);
% end;

Bt = Params.Bt;
Br = Params.Br;
% Gtrue = Params.Gtrue;

%% reconstruct the signal based on the quantized output
if bit == +inf
    r = r_bit;
else
    r = (- 2^( bit -1 ) + real(r_bit) + 1/2 ) * stepsize + ...
        1j * (- 2^( bit -1 ) + imag(r_bit) + 1/2 ) * stepsize;
    r = r + dither_mean + 1j* dither_mean;
end;

R = reshape(r, Nr, []);


% V_cyclic_0 = kron(eye(Nd), Bt') * T_cyclic;
% fast implementation of V_cyclic = kron(eye(Nd), Bt') * T_cyclic
T_temp = reshape(T_cyclic, Nt, Nd, Np);
V_temp = Bhfast(T_temp, Nt);
V_cyclic = reshape(V_temp, Nd*Nt, Np); 

Reps = 1;
tic;
% assume the gvar is 1/Nd
% Ghat_LMMSE = ctranspose( Br ) * R * V_cyclic' / ( (1-beta) * V_cyclic * V_cyclic' + ...
%     beta* Nd * SNR_linear * eye(Nt*Nd)+ Nd*eye(Nt*Nd) );

% use the estimated variance of g, i.e., Params.gvar
for ii=1:1:Reps
Ghat_LMMSE = ctranspose( Br ) * R * V_cyclic' / ( (1-beta) * V_cyclic * V_cyclic' + ...
    ( beta* SNR_linear * Nd + 1/Params.gvar ) * eye(Nt*Nd) );
end;
timeALMMSE = toc/Reps;

ghat_ALMMSE = Ghat_LMMSE(:);

Gtrue = Params.Gtrue;
% 10 * log10 ( (norm(ghat_LMMSE -  Gtrue(:)))^2/norm(Gtrue(:))^2 )

% tic;
% for ii=1:1:10
%     (1-beta) * V_cyclic * V_cyclic' + ...
%     ( beta* SNR_linear * Nd + 1/Params.gvar ) * eye(Nt*Nd);
% end;
% toc
% 
% tic;
% for ii=1:1:10
%     (1-beta) * V_cyclic * V_cyclic' + ...
%     diag(( beta* SNR_linear * Nd + 1/Params.gvar ) * ones(Nt*Nd, 1));
% end;
% toc
%Notes: The ALMMSE algorithm could be implemented by GAMP algorithm.
%However, GAMP algorithm may not be faster than this algorithm.

end

