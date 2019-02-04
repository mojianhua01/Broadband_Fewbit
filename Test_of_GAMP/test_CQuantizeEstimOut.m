close all;
clear vars;
dbstop if error;
M = 1e7; % # observations
tauw = 0.01; % observation noise variance
taup = 1; % prior variance on pre-quantized outputs

% generate realizations
%phat = linspace(-1,1,M).'; % prior mean (i.e., noisy estimates) of pre-quantized outputs
phat = abs(randn(M,2)) * [1; 1j] ; % prior mean (i.e., noisy estimates) of pre-quantized outputs
v = sqrt(taup/2)*randn(M,2) * [1 ; 1j];
z = phat + v; % pre-quantized outputs
w = sqrt(tauw/2)*randn(M,2) * [1; 1j]; % AWGN measurement noise

y = z + w;

Lloyd_stepsize_vector = [1.5956 0.9957 0.586 0.3352 0.1881 0.1041 0.0569 0.0308];
dither_mean = 0; % default value of dithering threshold

stepsize_scale_factor = 1;  % this stepsize setup works better for EMBGAMP and EMGMAMP algorithms
%                 stepsize_scale_factor = 1;  % this stepsize setup works better for BPDN algorithm

bit = 1;

stepsize = stepsize_scale_factor * rms([real(y); imag(y)])* Lloyd_stepsize_vector(bit);

% r_bit = min ( max( floor( (real(y)-dither_mean)/stepsize) + 2^(bit - 1), 0), 2^bit-1) ; %[0 2^bit-1]


r_bit_real = min ( max( floor( (real(y)-dither_mean)/stepsize) + 2^(bit - 1), 0), 2^bit-1) ; %[0 2^bit-1]
r_bit_imag = min ( max( floor( (imag(y)-dither_mean)/stepsize )+ 2^(bit - 1), 0), 2^bit-1) ; %[0 2^bit-1]
%
r_bit = r_bit_real + 1j * r_bit_imag;

%     r = sign(real(y)) .* ( min( ceil( abs(real(y)) /stepsize) , 2^(bit-1) ) - 1/2 ) * stepsize  + ...
%         1j* sign(imag(y)) .* ( min( ceil( abs(imag(y)) /stepsize) , 2^(bit-1) ) - 1/2 ) * stepsize;

% estimation under mismatch
mismatch_ = logspace(-0.1,0.1,5) % mismatch parameter
tauz_ = nan(size(mismatch_));
zerrmse_ = nan(size(mismatch_));
tauz_avg_ = nan(size(mismatch_));
for i=1:length(mismatch_)
    
    mismatch = mismatch_(i);
    
    if 1 % ideal case: {phat,taup} are correct
        phat_assumed = phat;
        taup_assumed = taup;
    else % this happens at first GAMP iteration (since xhat=0 => phat=0)
        phat_assumed = zeros(M,1);
        taup_assumed = mean(phat.^2) + taup;
    end
    
    EstimOut = CQuantizeEstimOut(r_bit, bit, stepsize, 0,tauw);
    [zhat,tauz] = EstimOut.estim(phat_assumed,mismatch*taup_assumed*ones(M,1));
    zerr = zhat-z; % estimation error
    zerrmse_(i) = mean(abs(zerr).^2); % MSE
    tauz_avg_(i) = mean(tauz); % predicted MSE
    
    figure(1); clf;
    hist( [real(zerr); imag(zerr)], 100)
    ylabel('histogram')
    xlabel('estimation error')
    grid on
    drawnow
    
end

figure(2); clf;
plot(mismatch_,zerrmse_,'o-', mismatch_,tauz_avg_,'x-')
legend('mse','predicted mse')
xlabel('scaling parameter (to test mismatch)')
grid on
