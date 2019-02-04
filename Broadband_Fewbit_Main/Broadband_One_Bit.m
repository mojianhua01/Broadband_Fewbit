clear;
close all;

% Model parameters
Nt = 8;
Nr = 8;
Nd = 16;
srrc = false; % use SRRC instead of RC time-domain pulse?
alf_rc = 0.5; % raised-cosine parameter [dflt=0.5]
Npre = min(2,Nd/2); % channel taps before first arrival [dflt=5]
Npst = min(2,Nd/2); % channel taps after last arrival [dflt=5]
Num_cluster = 2;   % number of path
plot_datacube = 1;
alpha = 10;
sparseRat = alpha * Num_cluster/Nd/Nt/Nr;

% Pilot length and SNR
N = 2^10;
SNR_vector = 10;


BGmean = 0;        	% Bernoulli-Gaussian active mean
BGvar = 1;         	% Bernoulli-Gaussian active variance

% create aperture/delay channel
H = zeros(Nr,Nt,Nd); % MIMO matrix for each delay
delay = zeros(Num_cluster,1);
gain = zeros(Num_cluster,1);
theta_t = zeros(Num_cluster,1);
theta_r = zeros(Num_cluster,1);
n = [0:Nd-1]+eps;               % baud-normalized sampling times
for l=1:Num_cluster % for each path...
    delay(l) = Npre+rand(1)*(Nd-Npre-Npst); % delay (between Npre and Nd-Npst)
    if srrc
        pulse = ( (1-alf_rc)*sinc((n-delay(l))*(1-alf_rc)) ...
            + cos(pi*(n-delay(l))*(1+alf_rc))*4*alf_rc/pi ...
            )./(1-(4*alf_rc*(n-delay(l))).^2);
    else % rc pulse
        pulse = cos(pi*alf_rc*(n-delay(l)))./(1-(2*alf_rc*(n-delay(l))).^2) ...
            .*sinc(n-delay(l)); % raised-cosine impulse response
    end
    pulse = pulse/norm(pulse);
    %plot(pulse,'.-'); title('pulse'); pause;
    
    theta_t(l) = rand(1); % transmit angle (between 0 and 1)
    theta_r(l) = rand(1); % receive angle (between 0 and 1)
    gain(l) = randn(1,2)*[1;1i]/sqrt(2*Num_cluster); % complex gain
    for d=1:Nd % for each delay ...
        H(:,:,d) =  H(:,:,d) + pulse(d)*gain(l)*sqrt(1/Nt/Nr) ...
            *exp(1i*2*pi*theta_r(l)*(0:Nr-1))' ...
            *exp(1i*2*pi*theta_t(l)*(0:Nt-1));
        Gtrue(:,:,d) = dftmtx(Nr)'/sqrt(Nr) * H(:,:,d) * dftmtx(Nt)/sqrt(Nt);
    end
end %l

gtrue = reshape(Gtrue, [], 1);

noise = (randn(Nr, Nt, N) + 1i * randn(Nr, Nt, N))/sqrt(2);


% Walsh_Hadamard_matrix = hadamard(N);
% Z_train = Walsh_Hadamard_matrix(randsample(1:N,Nt), 1:N);

% DFT_matrix = dftmtx(N);
% Z_train = DFT_matrix( randsample(1:N,Nt), 1:N );

% % Z_train = (randn(Nt, N) + 1j * randn(Nt, N));
Z_train = sign(randn(Nt,N))/sqrt(2)+1i*sign(randn(Nt,N))/sqrt(2);
% 
% % test
% Z_train = dftmtx(Nt)'*Z_train;

for index_SNR=1:1:length(SNR_vector)
    SNR = SNR_vector(index_SNR);
    SNR_linear = 10^(SNR/10);
    
    Z_train_scaled = Z_train./norm(Z_train, 'fro') *sqrt(N)*sqrt(SNR_linear);
    % scale Z_train such that sum((Z_train_scaled(:,ii)).^2) = SNR_linear for
    % ii=1,2,...
    
    % create circular-delay matrices: J{d} rotates down by (d-1) places
    J = cell(1,Nd-1);
    J{1} = speye(N);
    for d=1:Nd-1
        J{d+1} = [spalloc(N-d,d,0), speye(N-d); speye(d), spalloc(d,N-d,0)];
    end
    
    S = zeros(Nt*Nd, N);
    S(1: Nt,:) = Z_train_scaled;
    for d=2:1:Nd
        S( (d-1)*Nt+1: d*Nt, :) = Z_train_scaled * J{d};
    end;
    
    % A = kron(transpose(S), dftmtx(Nr)/sqrt(Nr));
    % z = A *gtrue;      % "True" transform vector
    
    z = reshape( dftmtx(Nr)/sqrt(Nr) * reshape(Gtrue, Nr, []) * S, [], 1);
    
    wvar = 1;
    
    % Compute training class labels based on synthetic "true" hyperplane
    %     y_train = double(z > 0);
    %         r = z + noise/sqrt(SNR_linear);
    y = z + sqrt(1/2)*wvar*randn(Nr*N,2)*[1;1i];
    %         r = z + sqrt(1/2)*randn(400,1)/sqrt(SNR_linear);
    r_train = sign(real(y)) + 1j* sign(imag(y));
    
    
    % Decide on MAP or MMSE GAMP
    map = false;
    
    % Create an input estimation class corresponding to a Gaussian vector
    % inputEst0 = CAwgnEstimIn(BGmean, BGvar, map);
    % inputEst = SparseScaEstim(inputEst0, sparseRat);
    
    % Auto tune the mean, variance and sparseRat.
    % This is extremely important for the algorithm to work!!!
    inputEst = CAwgnEstimIn(BGmean, BGvar, map, 'autoTune', true);
    inputEst = SparseScaEstim(inputEst, sparseRat, 0, 'autoTune', true);
    
    % Create an output estimation class corresponding to the Gaussian noise.
    % Note that the observation vector is passed to the class constructor, not
    % the gampEst function.
    outputEst = CProbitEstimOut(r_train, 0,wvar, map);
    % outputEst = CQuantizeEstimOut(y_train, 1, 1, 0, wvar, map);
    
    % Set the default options
    opt = GampOpt();
    opt.nit = 20;
    opt.tol = 1e-4;
    opt.uniformVariance=0;
    opt.pvarMin=0;
    opt.xvarMin=0;
    opt.adaptStep =  true;
    opt.adaptStepBethe= false;
    opt.legacyOut=false;
    opt.xvarMin = 1e-4;     % At very high SNR, use very small pvarMin!
    % Demonstrate automatic selection of xvar0
    % if 1
    %   opt.xvar0auto = true;
    %   opt.xhat0 = x + 1*(randn(nx,2)*[1;1j])*norm(x)/sqrt(nx);
    % end;
    
    hA = @(in) Afast_Digital(in,S,Nr,1);
    hAt = @(in) Afast_Digital(in,S,Nr,0);
    lA = FxnhandleLinTrans(Nr*N, Nt*Nr*Nd, hA, hAt); % linear transform object
    
    % Run the GAMP algorithm
    tic
    [estFin,optFin,estHist] = gampEst(inputEst, outputEst, lA, opt);
    xhat = estFin.xhat;
    timeGAMP = toc
    
    GAMP_All(index_SNR) = norm(xhat(:, end) - gtrue)^2;
end;

GAMP_NMSE = GAMP_All./norm(gtrue)^2;
figure,
plot(SNR_vector, 10*log10(GAMP_NMSE), 'r');
xlabel('SNR', 'fontsize',14)
ylabel('Normalized MSE', 'fontsize',14)

% Now perform the exact LMMSE solution
xmean0 = zeros(Nt*Nr*Nd,1);
xvar0 = BGvar* ones(Nt*Nr*Nd,1);
% tic
% xhatLMMSE = xmean0 + xvar0.*(A'*((A*diag(xvar0)*A'+diag(wvar))\(r-A*xmean0)));
% timeLMMSE = toc

% Plot the results
figure;
[xsort,I] = sort(real(gtrue));
plot(xsort, xsort,'-', xsort,[real(xhat(I))], '.');
set(gca,'FontSize',16);
grid on;
legend('True', 'GAMP estimate');
xlabel('True value of real-x');
ylabel('Estimate of real-x');

figure;
subplot(311)
for ii = 2:1:size(estHist.xhat, 2)
    MSE(ii) = 10*log10(sum(abs( estHist.xhat(:,ii)/sqrt(sum(abs(estHist.xhat(:,ii)).^2)) - ...
        gtrue/norm(gtrue)^2 ).^2));
end;
plot(MSE);
ylabel('NMSE [dB]')
grid on
subplot(312)
plot(estHist.step)
ylabel('step')
grid on
subplot(313)
plot(estHist.val)
ylabel('val')
xlabel('iteration')
grid on

figure,
subplot(121)
gtrue2 = reshape(permute(reshape(gtrue,Nr,Nt,Nd),[1,3,2]),Nd*Nr,Nt);
bar3(abs(gtrue2));
title('True value of the channel', 'Fontsize', 13)
xlabel('Tx angle')
ylabel('Rx angle & delay')

subplot(122)
ghat2 = reshape(permute(reshape(estHist.xhat(:,end),Nr,Nt,Nd),[1,3,2]),Nd*Nr,Nt);
bar3(abs(ghat2));
title('Estimate of the channel', 'Fontsize', 13)
xlabel('Tx angle')
ylabel('Rx angle & delay')

% Display the MSE
mseGAMP = 20*log10( norm(gtrue-xhat)/norm(gtrue));
% mseLMMSE = 20*log10( norm(x-xhatLMMSE)/norm(x));
fprintf(1,'GAMP:  MSE = %5.1f dB\n', mseGAMP);
% fprintf(1,'LMMSE: MSE = %5.1f dB\n', mseLMMSE);

% figure,
% bar3(abs(gtrue2-ghat2));
% title('Estimation Error', 'Fontsize', 13)
% xlabel('Tx angle')
% ylabel('Rx angle & delay')

% compute the coefficient between the real and imaginary parts
correlation = sum(real(gtrue).* imag(gtrue))/norm(real(gtrue))/norm(imag(gtrue))