function [ ghat_Lasso ] = func_Broadband_Few_Bit_UPA_Lasso( T_cyclic, r_bit, Params)
% Estimate the channel given the training signals and quantization outputs
% Inputs:
%   T_cyclic:   training signals
%   r_bit:      Quantization resolution
%   Params:     Parameters
% Output:
%   ghat:       estimate of the channel in the angular domain

% Channel model
Nt = Params.Nt;
Nta = Params.Nta;
Nte = Params.Nte;
Nr = Params.Nr;
Nra = Params.Nra;
Nre = Params.Nre;
Nd = Params.Nd;
Np = Params.Np;

bit = Params.bit;
stepsize = Params.stepsize;
% SNR_linear = Params.SNR_linear;
dither_mean = Params.dither_mean;
% wvar = Params.wvar;

Bt = Params.Bt;

% V_cyclic_0 = kron(eye(Nd), Bt') * T_cyclic;

% % fast implementation of V_cyclic = kron(eye(Nd), Bt') * T_cyclic
% T_temp = reshape(T_cyclic, Nt, Nd, Np);
% V_temp = Bhfast(T_temp, Nt);
% V_cyclic_1 = reshape(V_temp, Nd*Nt, Np);

V_cyclic = NaN(Nr*Nd, Np);
% V_train = Bt'* T_cyclic(1:Nt,:);
V_train = Bhfast(T_cyclic(1:Nt,:), Nt);
for ii=0:1:Nd-1
    V_cyclic( ii*Nt+1: (ii+1)*Nt, :) = circshift(V_train,ii,2);
end

% Another implementation of V_cyclic which has similar efficiency.
% T_temp = reshape(T_cyclic, Nt, Nd, Np);
% V_temp = Bhfast(T_temp, Nt);
% V_cyclic_2 = reshape(V_temp, Nd*Nt, Np); 

%% MATLAB BPDN
% X = ( 1 - beta ) * kron( S.', 1/sqrt(Nr) * kron(dftmtx(Nra), dftmtx(Nre)) );
% Y = r;
% Ghat_BPDN = lasso(X, Y);

% A = ( 1 - beta ) * kron( S.', 1/sqrt(Nr) * kron( dftmtx(Nra), dftmtx(Nre) ) );
% sigma = 1.0 * sqrt( (1-beta) * wvar + beta * (1-beta) * SNR_linear );
% opts = spgSetParms('verbosity',0);
% Ghat_BPDN = spg_bpdn(A, r, sigma, opts);

% Lloyd_stepsize_vector = [1.5956 0.9957 0.586 0.3352 0.1881 0.1041 0.0569 0.0308];
% beta_vector = [0.3634 0.1188 0.03744 0.01154];


%% Compute mu (mu = 1 - beta) and sigma by simulations.  
% only consider the real part, so the noise variance is 1/2
if bit == +inf
    mu = 1;
    sigma = sqrt(1/2);
else
    % r = Q(z + w) = mu * (z + w) + noise = mu * z + ( mu * w + noise )
%     stepsize_scale_factor = Params.stepsize_scale_factor;
    y_temp = randn(10000,1)*sqrt(Params.received_power_linear);
    stepsize_temp =  stepsize;
    r_temp = sign( y_temp ) .* ( min( ceil( abs( y_temp )/stepsize_temp) , 2^(bit-1) ) - 1/2 ) * stepsize_temp;
    mu = mean(y_temp.*r_temp)/rms(y_temp)^2;
    sigma = sqrt ( rms((r_temp - mu.*y_temp))^2 + mu^2 * 1/2 );
%     corr = mean(z_temp.*r_temp)/rms(z_temp)/rms(r_temp)
end;

%% reconstruct the signal based on the quantized output
if bit == +inf
    r = r_bit;
else
    r = (- 2^( bit -1 ) + real(r_bit) + 1/2 ) * stepsize + ...
        1j * (- 2^( bit -1 ) + imag(r_bit) + 1/2 ) * stepsize;
    r = r + dither_mean + 1j* dither_mean;
end;

%% spgl1 BPDN
% tic;
% % A = ( 1 - beta ) * kron( S.', 1/sqrt(Nr) * kron( dftmtx(Nra), dftmtx(Nre) ) );
% opA = @(in, mode) (1-beta)*Afast_Digital_UPA(in, T_cyclic, V_cyclic, Params, mode);
% sigma =  1.1* sqrt( length(r) * ((1-beta) * wvar + beta * (1-beta) * SNR_linear ));
% opts = spgSetParms('verbosity',0);
% ghat_BPDN = spg_bpdn(opA, r, sigma, opts);
% 
% %Ghat_BPDN = reshape(Ghat_BPDN, Nr, Nt, Nd);
% 
% timeBPDN = toc;

% tic;
% A = ( 1 - beta ) * kron( S.', 1/sqrt(Nr) * kron( dftmtx(Nra), dftmtx(Nre) ) );
opA = @(in, mode) mu * Afast_Digital_UPA(in, T_cyclic, V_cyclic, Params, mode);

coeff = 2*6/pi^2; % for H_UPA_16_4_16_2clusters
% coeff = 4*6/pi^2; % for H_UPA_64_64_16_4clusters
tau = Params.estimated_norm * log(Nt*Nr*Nd) * coeff;

Gnorm_1 = norm(Params.Gtrue(:),1);

opts = spgSetParms('verbosity',1);
ghat_Lasso = spg_lasso(opA, r, tau, opts);

%Ghat_Lasso = reshape(Ghat_BPDN, Nr, Nt, Nd);

% timeLasso = toc;

end

