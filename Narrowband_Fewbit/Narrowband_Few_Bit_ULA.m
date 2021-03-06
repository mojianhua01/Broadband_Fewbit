clear;
clc;
close all;
% Model parameters
Nt = 32;
Nr = 32;
Num_cluster = 2;
Num_ray = 10;
Angle_spread = pi/24; % pi/24;

[~, x0] = func_mmWave_Channel_ULA(Nt, Nr, Num_cluster, Num_ray, Angle_spread);

x = reshape(x0,[],1);

L = 2^12;

bit = 2;

N = Nt*Nr;           % Number of features (dimension of x)
sparseRat = 10 * Num_cluster/Nt/Nr;

BGmean = 0;        	% Bernoulli-Gaussian active mean
BGvar = 1;         	% Bernoulli-Gaussian active variance

SNR_vector = 20;
for index_SNR=1:1:length(SNR_vector)
    SNR = SNR_vector(index_SNR);
    SNR_linear = 10^(SNR/10);
    %% Walsh Hadamard Training Signals
    % Walsh_Hadamard_matrix = hadamard(L);
    % Z_train = Walsh_Hadamard_matrix(randsample(1:L,Nt), 1:L);
    
    %% Gaussian Random Training Signals
    % Z_train = (randn(Nt, L) + 1j* randn(Nt, L));
    
    %% 1/-1 Training Signals
    Z_train = sign(randn(Nt,L))/sqrt(2)+1i*sign(randn(Nt,L))/sqrt(2);
    
    %% scale Z_train such that sum((Z_train_scaled(:,ii)).^2) = SNR_linear for ii=1,2,...L
    Z_train_scaled = Z_train./sqrt(sum(sum(abs(Z_train).^2)))*sqrt(L)*sqrt(SNR_linear);
    
    
    % A = kron(transpose(Z_train_scaled), dftmtx(Nr)/sqrt(Nr));  % A_train_cv is complex-valued
    % z = A *x;      % "True" transform vector
    
    global z;
    z = reshape( dftmtx(Nr)/sqrt(Nr) * x0 * Z_train_scaled, [], 1);
    
    wvar = 1;
    
    % Compute training class labels based on synthetic "true" hyperplane
    %     y_train = double(z > 0);
    %         r = z + noise/sqrt(SNR_linear);
    y = z + sqrt(wvar/2)*randn(Nr*L,2)*[1;1i];
    
    if bit == 1;
        stepsize = max(abs([real(y); imag(y)])) * 1.01;  % if ADC has only one bit, set the stepsize correctly
    else
        stepsize = 1.01*( max(abs([real(y); imag(y)])))/(2^(bit-1));
    end;
    
    r_bit_real = min ( max( floor( real(y)/stepsize) + 2^(bit - 1), 0), 2^bit-1) ;
    r_bit_imag = min( max( floor( imag(y)/stepsize )+ 2^(bit - 1), 0), 2^bit-1);
    
    r_bit = r_bit_real + 1j * r_bit_imag;
    
    % probit_var = 1/SNR_linear;
    
    hA = @(in) Afast_Digital(in,Z_train_scaled,Nr,1);
    hAt = @(in) Afast_Digital(in,Z_train_scaled,Nr,0);
    lA = FxnhandleLinTrans(Nr*L,Nt*Nr,hA,hAt); % linear transform object
    
    % Decide on MAP or MMSE GAMP
    map = false;
    
    % Create an input estimation class corresponding to a Gaussian vector
    % inputEst0 = CAwgnEstimIn(BGmean, BGvar, map);
    % inputEst = SparseScaEstim(inputEst0, sparseRat);
    
    % Auto tune the mean, variance and sparseRat.
    % This is extremely important for the algorithm to work!!!
    inputEst0 = CAwgnEstimIn(BGmean, BGvar, map, 'autoTune', true);
    inputEst = SparseScaEstim(inputEst0, sparseRat, 0, 'autoTune', true);
    
    % Create an output estimation class corresponding to the Gaussian noise.
    % Note that the observation vector is passed to the class constructor, not
    % the gampEst function.
    outputEst = CQuantizeEstimOut(r_bit, bit, stepsize, 0, wvar, map);
    
    % Set the default options
    opt = GampOpt();
    opt.nit = 15;
    opt.tol = 1e-4;
    opt.uniformVariance=0;
    opt.pvarMin=1e-5;
    opt.xvarMin=1e-5;
    opt.adaptStep=true;
    opt.adaptStepBethe=false;
    opt.legacyOut=false;
    opt.verbose=1;
    opt.step=1;
    

    
    % Demonstrate automatic selection of xvar0
    % if 1
    %   opt.xvar0auto = true;
    %   opt.xhat0 = x + 1*(randn(length(x),2)*[1;1j])*norm(x)/sqrt(length(x));
    % end;
    
    % Run the GAMP algorithm
    tic
    [estFin,optFin,estHist] = gampEst(inputEst, outputEst, lA, opt);  % fast operator
    % [estFin,optFin,estHist] = gampEst(inputEst, outputEst, A, opt);
    xhat = estFin.xhat;
    timeGAMP = toc
    
    GAMP_All(index_SNR) = norm(xhat(:, end) - x)^2;
end;

GAMP_NMSE = GAMP_All./norm(x)^2;
figure,
plot(SNR_vector, 10*log10(GAMP_NMSE), 'r');
xlabel('SNR', 'fontsize',14)
ylabel('Normalized MSE', 'fontsize',14)

% Now perform the exact LMMSE solution
xmean0 = zeros(Nt*Nr,1);
xvar0 = BGvar* ones(Nt*Nr,1);
% tic
% xhatLMMSE = xmean0 + xvar0.*(A'*((A*diag(xvar0)*A'+diag(wvar))\(r-A*xmean0)));
% timeLMMSE = toc;

%%convex method proposed by Junil Choi. It is possible when A_train
%%has small dimension, but will not work when A_train has large
%%dimensions.
% signed_A = diag(y_train)*A;
% signed_A = [real(signed_A), -imag(signed_A) ; imag(signed_A) real(signed_A)];
% cvx_begin
% variable x(2*Nr*Nt);
% maximize sum(log_normcdf(signed_A*x));
% subject to
% norm(x) <= 2*Nt*Nr;
% cvx_end

% % Plot the results
% figure;
% [xsort,I] = sort(real(x));
% % plot(xsort, xsort,'-', xsort,[real(xhat(I)) real(xhatLMMSE(I))], '.');
% plot(xsort, xsort,'-', xsort,[real(xhat(I))], '.');
% set(gca,'FontSize',16);
% grid on;
% % legend('True', 'GAMP estimate', 'LMMSE estimate');
% legend('True', 'GAMP estimate');
% xlabel('True value of real-x');
% ylabel('Estimate of real-x');

% figure;
% subplot(311)
% plot(10*log10(sum(abs( estHist.xhat - x*ones(1,size(estHist.xhat,2)) ).^2,1)/norm(x)^2))
% ylabel('NMSE [dB]')
% grid on
% subplot(312)
% plot(estHist.step)
% ylabel('step')
% grid on
% subplot(313)
% plot(estHist.val)
% ylabel('val')
% xlabel('iteration')
% grid on

% Display the MSE
mseGAMP = 20*log10( norm(x-xhat)/norm(x));
% mseLMMSE = 20*log10( norm(x-xhatLMMSE)/norm(x));
fprintf(1,'GAMP:  MSE = %5.1f dB\n', mseGAMP);
fprintf(1,'GAMP:  Time = %5.1f sec\n', timeGAMP);
% fprintf(1,'LMMSE: MSE = %5.1f dB\n', mseLMMSE);

% figure,
% for ii = 2:1:size(estHist.xhat, 2)
%     MSE(ii) = 10*log10(sum(abs( estHist.xhat(:,ii)/sqrt(sum(abs(estHist.xhat(:,ii)).^2)) - ...
%         x/norm(x)^2 ).^2));
% end;
% plot(MSE);
% ylabel('NMSE [dB]')

figure,
subplot(121)
bar3(abs(x0))
title('True value of the channel', 'Fontsize', 14)
xlabel('Tx angle','Fontsize', 12)
ylabel('Rx angle','Fontsize', 12)
subplot(122)
Estimation = reshape(estHist.xhat(:,end), Nr, Nt);
bar3(abs(Estimation))
title('Estimate of the channel', 'Fontsize', 14)
xlabel('Tx angle','Fontsize', 12)
ylabel('Rx angle','Fontsize', 12)