function [ ghat_LS ] = func_Broadband_Few_Bit_UPA_LS( T_cyclic, r_bit, Params)
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
SNR_linear = Params.SNR_linear;
dither_mean = Params.dither_mean;
wvar = Params.wvar;

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
for ii=1:1:Reps
% Ghat_LS = ctranspose( 1/sqrt(Nr) * kron(dftmtx(Nra), dftmtx(Nre)) ) * R * V' * ( V * V')^(-1);
Ghat_LS = ctranspose( Br ) * R * V_cyclic' / ( V_cyclic * V_cyclic');
end;
timeLS = toc/Reps;
ghat_LS = Ghat_LS(:);

Gtrue = Params.Gtrue;
% 10*log10 ( (norm(ghat_LS -  Gtrue(:)))^2/norm(Gtrue(:))^2 );

end

