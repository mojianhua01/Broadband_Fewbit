clear;
% Model parameters
Nt = 64;
Nr = 16;
L = 256;
alpha = 10;
Num_cluster = 2;

% load Data.mat
load H_64_16_2paths.mat;
index_H = 1;
x0 = Hv_all(1:Nr,1:Nt, index_H);
%     x0 = H_all(1:Nr, 1:Nt, index_H);
%     x0 = x0.*(rand(Nr, Nt) < sparseRat);
x = reshape(x0,[],1);

noise_cv = reshape(Noise(1:Nr, 1:L), [], 1);

% [H_scaled, H_v] = func_mmWave_Channel_UPA(BS_ant_W, BS_ant_H, MS_ant_W,MS_ant_H, Num_cluster, Num_ray, Angle_spread)

N = Nt*Nr;           % Number of features (dimension of x)
sparseRat = alpha * Num_cluster/Nt/Nr;

BGmean = 0;        	% Bernoulli-Gaussian active mean
BGvar = 1;         	% Bernoulli-Gaussian active variance


SNR = 10;
SNR_linear = 10^(SNR/10);
% Walsh_Hadamard_matrix = hadamard(L);
% Z_train = Walsh_Hadamard_matrix(randsample(1:L,Nt), 1:L);

Z_train = (randn(Nt, L) + 1j* randn(Nt, L));
Z_train_scaled = Z_train./sqrt(sum(sum(abs(Z_train).^2)))*sqrt(L);
% scale Z_train such that sum((Z_train_scaled(:,ii)).^2) = 1 for
% ii=1,2,...L

A = kron(transpose(Z_train_scaled), dftmtx(Nr)/sqrt(Nr));  % A_train_cv is complex-valued
z = A *x;      % "True" transform vector

wvar = 1/SNR_linear;

% Compute training class labels based on synthetic "true" hyperplane
%     y_train = double(z > 0);
%         r = z + noise/sqrt(SNR_linear);
r = z + sqrt(wvar/2)*randn(Nr*L,2)*[1;1i];
%         r = z + sqrt(1/2)*randn(400,1)/sqrt(SNR_linear);
y_train = sign(real(r)) + 1j* sign(imag(r));

probit_var = 1/SNR_linear;

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
outputEst = CProbitEstimOut(y_train, 0, wvar, map);

% Set the default options
opt = GampOpt();
opt.nit = 1000;
opt.tol = 1e-4;
opt.uniformVariance=0;
opt.pvarMin=0;
opt.xvarMin=0;
opt.adaptStep=true;
opt.adaptStepBethe=true;
opt.legacyOut=false;

% Demonstrate automatic selection of xvar0
% if 1
%   opt.xvar0auto = true;
%   opt.xhat0 = x + 1*(randn(length(x),2)*[1;1j])*norm(x)/sqrt(length(x));
% end;

% Run the GAMP algorithm
tic
[estFin,optFin,estHist] = gampEst(inputEst, outputEst, A, opt);
xhat = estFin.xhat;
timeGAMP = toc;

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

% Plot the results
figure;
[xsort,I] = sort(real(x));
% plot(xsort, xsort,'-', xsort,[real(xhat(I)) real(xhatLMMSE(I))], '.');
plot(xsort, xsort,'-', xsort,[real(xhat(I))], '.');
set(gca,'FontSize',16);
grid on;
% legend('True', 'GAMP estimate', 'LMMSE estimate');
legend('True', 'GAMP estimate');
xlabel('True value of real-x');
ylabel('Estimate of real-x');

figure;
subplot(311)
plot(10*log10(sum(abs( estHist.xhat - x*ones(1,size(estHist.xhat,2)) ).^2,1)/norm(x)^2))
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

% Display the MSE
mseGAMP = 20*log10( norm(x-xhat)/norm(x));
% mseLMMSE = 20*log10( norm(x-xhatLMMSE)/norm(x));
fprintf(1,'GAMP:  MSE = %5.1f dB\n', mseGAMP);
fprintf(1,'GAMP:  Time = %5.1f sec\n', timeGAMP);
% fprintf(1,'LMMSE: MSE = %5.1f dB\n', mseLMMSE);

subplot(121)
bar3(abs(x0))
subplot(122)
Estimation = reshape(estHist.xhat(:,end), Nr, Nt);
bar3(abs(Estimation))