function [ Ghat_LMMSE ] = func_Broadband_Few_Bit_ULA_LMMSE( Gtrue, SNR, S_train, bit)
% Need revision!!!
% Estimate the channel given the SNR and training signals
% Inputs:
%   S_train:    Nr by N training signals
%   Gtrue:      true value of the channel in the angular domain
%   Bit:        Quantization resolution
%   SNR:        Transmission power in dB
%
% Output:
%   Ghat:       estimate of the channel in the angular domain

% Channel model
Nt = size(Gtrue, 2);
Nr = size(Gtrue, 1);
Nd = size(Gtrue, 3);
N = size(S_train, 2);

SNR_linear = 10^(SNR/10);

S_train_scaled = S_train./(norm(S_train, 'fro'))*sqrt(N)*sqrt(SNR_linear);
% scale S_train such that sum((S_train_scaled(:,ii)).^2) = SNR_linear for
% ii=1,2,...

% create circular-delay matrices: J{d} rotates down by (d-1) places
J = cell(1,Nd-1);
J{1} = speye(N);
for d=1:Nd-1
    J{d+1} = [spalloc(N-d,d,0), speye(N-d); speye(d), spalloc(d,N-d,0)];
end

S = zeros(Nt*Nd, N);
S(1: Nt,:) = S_train_scaled;
for d=2:1:Nd
    S( (d-1)*Nt+1: d*Nt, :) = S_train_scaled * J{d};
end;

% A = kron(transpose(S), dftmtx(Nr)/sqrt(Nr));
% z = A *gtrue;      % "True" transform vector

z = reshape( dftmtx(Nr)/sqrt(Nr) * reshape(Gtrue, Nr, []) * S, [], 1);

wvar = 1;

% Compute training class labels based on synthetic "true" hyperplane
y = z + sqrt(1/2)*wvar*randn(Nr*N,2)*[1;1i];

stepsize_vector = [1.5956 0.9957 0.586 0.3352 0.1881 0.1041 0.0569 0.0308];
beta_vector = [0.3634 0.1188 0.03744 0.01154];
if bit == +inf
    r = y;
else
    %     if bit == 1;
    %         stepsize = 1.01 * max(abs([real(y); imag(y)]));  % if ADC has only one bit, set the stepsize correctly
    %     else
    %         stepsize = ( 1.01 * max(abs([real(y); imag(y)])))/(2^(bit-1));
    %     end;
    
    if bit == 1;
        stepsize = rms([real(y); imag(y)])* stepsize_vector(bit);
        %     stepsize = max(abs([real(y); imag(y)])) * 1.01;  % if ADC has only one bit, set the stepsize correctly
    else
        stepsize = rms([real(y); imag(y)])* stepsize_vector(bit);
    end;
    
    r_bit_real = min ( max( floor( real(y)/stepsize) + 2^(bit - 1), 0), 2^bit-1) ;
    r_bit_imag = min( max( floor( imag(y)/stepsize )+ 2^(bit - 1), 0), 2^bit-1);
    
    r_bit = r_bit_real + 1j * r_bit_imag;
    
    r = (- 2^( bit -1 ) + r_bit_real + 1/2 ) * stepsize + ...
        1j * (- 2^( bit -1 ) + r_bit_imag + 1/2 ) * stepsize;
end;

R = reshape(r, Nr, []);
if bit == +inf
    beta = 0;
else
    beta = beta_vector(bit);
end;

tic;
% Ghat_LMMSE = ctranspose(dftmtx(Nr)/sqrt(Nr)) * R * S' * ( (1-beta) * S * S' + ...
%     beta* Nd * SNR_linear * eye(Nt*Nd)+ Nd*eye(Nt*Nd) )^(-1);
Ghat_LMMSE = ctranspose(dftmtx(Nr)/sqrt(Nr)) * R * S' / ( (1-beta) * S * S' + ...
    beta* Nd * SNR_linear * eye(Nt*Nd)+ Nd*eye(Nt*Nd) );
Ghat_LMMSE = reshape(Ghat_LMMSE, Nr, Nt, Nd);
timeLMMSE = toc;

end

