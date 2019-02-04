function [ Ghat_QIST ] = func_Broadband_Few_Bit_UPA_QIST( T_cyclic, r_bit, Params)

% Have not finished

% Estimate the channel given the SNR and training signals
% Inputs:
%   T_train:    Nr by N training signals
%   Gtrue:      true value of the channel in the angular domain
%   Bit:        Quantization resolution
%   SNR:        Transmission power in dB
%
% Output:
%   Ghat:       estimate of the channel in the angular domain

% Channel model
Nt = Params.Nt;
Nta = Params.Nta;
Nte = Params.Nte;
Nr = Params.Nr;
Nra = Params.Nra;
Nre = Params.Nre;
Nd = Params.Nd;
Np = Params.Np;

SNR_linear = Params.SNR_linear;
bit = Params.bit;
stepsize = Params.stepsize;
% SNR_linear = Params.SNR_linear;
dither_mean = Params.dither_mean;
% wvar = Params.wvar;

Bt = Params.Bt;

V_cyclic = NaN(Nr*Nd, Np);
% V_train = Bt'* T_cyclic(1:Nt,:);
V_train = Bhfast(T_cyclic(1:Nt,:), Nt);
for ii=0:1:Nd-1
    V_cyclic( ii*Nt+1: (ii+1)*Nt, :) = circshift(V_train,ii,2);
end

stepsize_vector = [1.5956 0.9957 0.586 0.3352 0.1881 0.1041 0.0569 0.0308];
beta_vector = [0.3634 0.1188 0.03744 0.01154];


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

opA = @(in, mode) Afast_Digital_UPA(in, T_cyclic, V_cyclic, Params, mode);
quantize = @(y, stepsize) sign(real(y)) .* ( min( ceil( abs(real(y)) /stepsize) , 2^(bit-1) ) - 1/2 ) * stepsize  + ...
    1j* sign(imag(y)) .* ( min( ceil( abs(imag(y)) /stepsize) , 2^(bit-1) ) - 1/2 ) * stepsize;
ST = @(a,tau) sign(real(a)).*(max(abs(real(a))-tau, 0)) + 1j* sign(imag(a)).*(max(abs(imag(a))-tau, 0));
% tau = sigma/(mu^2 * SNR_linear * Nd); % tau is too large!!!
tau = 0.1* sqrt(2) * sigma/(mu * sqrt(SNR_linear * Np/Nt));
maxiter = 100;
stopcrit = 10^(-2);
xn(:,1) = zeros(Nt*Nr*Nd,1);
alpha = 0.1/(mu^2 * SNR_linear * Np/Nt);
%% Quantized iterative soft thresholding
tic;
for n = 1:maxiter
    
%     oxn = xn;
    
    % main QIST iteration
    
    Axn = opA(xn(:,n), 1);
    if bit == +inf
        QAxn = Axn;
    else
%         stepsize = stepsize_vector(bit) * rms(Axn);
        QAxn = quantize(Axn, stepsize);
    end;
    
    a = xn(:,n) + alpha * opA(r - QAxn, 2);
%     a = xn + alpha * opA(r - Axn, 2);
%     a = xn + alpha * A' * (r - Q(A*xn)); 
    xn(:,n+1) = ST(a,tau);
    
    if norm(xn(:,n+1)-xn(:,n))/norm(xn(:,n)) < stopcrit
        break;
    end
end
timeQIST = toc;

if 0
    Gtrue = Params.Gtrue;
    figure,
    plot(sum(abs(bsxfun(@minus, xn, Gtrue(:))).^2)./norm(Gtrue(:))^2)
    sum(abs(bsxfun(@minus, xn, Gtrue(:))).^2)/norm(Gtrue(:))^2;
    xlabel('Iteration')
    ylabel('Normalized MSE')
    title('QIST')
    grid on;
end;

Ghat_QIST = xn(:,end);

end

