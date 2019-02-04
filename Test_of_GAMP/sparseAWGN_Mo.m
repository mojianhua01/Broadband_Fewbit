% sparseAWGN:  Example of estimating a sparse vector with Gaussian noise.
%
% In this problem, x is a Bernoulli-Gaussian random vector, that we want to
% estimate from measurements of the form
%
%   y = A*x + w,
%
% where A is a random matrix and w is Gaussian noise.  This is a classical
% compressed sensing problem of estimating a sparse vector for random
% linear measurements.

clc;
clear;
% close all;
Nt = 10;
Nr = 10;
L = 256;
sparseRat = 0.05;
load Data.mat
load H_v.mat

% Add to path
addpath('../main/');

SNR_vector = -10:2:20;
for index_H = 1:1:1
    index_H
    H = H_all(1:Nr, 1:Nt, index_H);
    H_Sparse = H.*(Sparse_num(1:Nr, 1:Nt) < sparseRat);
    %     H_Sparse = H_v;
    h_cv = reshape(H_Sparse, [], 1);
    h = [real(h_cv); imag(h_cv)];
    noise_cv = reshape(Noise(1:Nr, 1:L), [], 1);
    noise = [real(noise_cv); imag(noise_cv)];

    for index_SNR=1:1:length(SNR_vector)
        SNR = SNR_vector(index_SNR);
        SNR_linear = 10^(SNR/10);

        %         Walsh_Hadamard_matrix = hadamard(L);
        %         Z_train_scaled = sqrt(SNR_linear/Nt)*Walsh_Hadamard_matrix(1:Nt, 1:L);

        Z_train = (randn(Nt, L) + 1j* randn(Nt, L));
        Z_train_scaled = Z_train./sqrt(sum(sum(abs(Z_train).^2)))*sqrt(L)*sqrt(SNR_linear);

        A_train_cv = kron(transpose(Z_train_scaled),  dftmtx(Nr)/sqrt(Nr));
        A_train = [real(A_train_cv) -imag(A_train_cv) ;
            imag(A_train_cv) real(A_train_cv)];

        r = A_train * h + noise;
        y = sign(r);
        %         lambda_vector = [0.1:0.2:0.5];
        %         for ii = 1:1:length(lambda_vector)
        %             lambda = lambda_vector(ii);
        %             tic;
        %             hhat = Lasso( A_train, r, lambda );
        %             toc;
        %             Lasso_MSE_All(index_H, index_SNR , ii) = (norm(hhat-h))^2/(norm(h)^2);
        %         end;

        % ML algorithm
        % hhat = (A_train'*A_train)^(-1)*A_train'*r;
        %         hhat = conjgrad(A_train'*A_train, A_train'*r, zeros(size(A_train,2),1));
        %         ML_MSE_All(index_H, index_SNR) = (norm(hhat(:,end)-h))^2/(norm(h)^2);

        % OMP algorithm
        %         hhat = OMP(A_train, r, 2*Nr*Nt*sparseRat);
        %         OMP_MSE_All(index_H, index_SNR) = (norm(hhat(:,end)-h))^2/(norm(h)^2);

        %LMMSE has similar performance as ML
        %         hhat = (A_train'*sparseRat*A_train + 1/2* eye())^(-1)*A_train'*sparseRat*r;
        %         LMMSE_MSE_All(index_H, index_SNR) = (norm(hhat-h))^2/(norm(h)^2);

        % Generate input estimation class
        % First, create an estimator for a Gaussian random variable (with no
        % sparsity)
        xmean0 = 0;
        xvar0 = 1;
        inputEst0 = AwgnEstimIn(xmean0, xvar0);

        % Then, create an input estimator from inputEst0 corresponding to a random
        % variable x that is zero with probability 1-sparseRat and has the
        % distribution of x in inputEst0 with probability sparseRat.
        inputEst = SparseScaEstim( inputEst0, sparseRat );

        % Output estimation class:  Use the
        outputEst = AwgnEstimOut(r, 1/2);

        % Run the GAMP algorithm
        % Set the default options
        opt = GampOpt();
        opt.nit = 10000;
        opt.tol = max(min(1e-3, 10^(-SNR/10)),1e-15);
        opt.uniformVariance=0;
        opt.pvarMin=0;
        opt.xvarMin=0;
        opt.adaptStep=true;
        opt.adaptStepBethe=true;
        opt.legacyOut=false;
        %         [xhat, xvar] = gampEst(inputEst, outputEst, A_train, GAMPopt);
        %         GAMP_no_quantize_MSE_All(index_H, index_SNR) = (norm(xhat-h))^2/(norm(h)^2);

        % Create the EstimOut class object for the probit channel (MMSE)
        probit_var = 1/2;
        maxSumVal = false;
        outEst = ProbitEstimOut(y, 0, probit_var, maxSumVal);
        [estFin,optFin,estHist] = gampEst(inputEst, outEst, A_train, opt);
        xhat = estFin.xhat;
        GAMP_SE_All(index_H, index_SNR) = norm(xhat-h)^2;
        h_power_All(index_H, index_SNR) = norm(h)^2;
        %         GAMP_MSE_All(index_H, index_SNR) = (norm(xhat-h))^2/(norm(h)^2);
    end;
end;
% figure,
% ML_MSE = mean(ML_MSE_All,1);
% plot(SNR_vector, 10*log10(ML_MSE));
% hold on;

% LMMSE_MSE = mean(LMMSE_MSE_All)
% plot(SNR_vector, 10*log10(LMMSE_MSE), 'r');

% OMP_MSE = mean(OMP_MSE_All);
% plot(SNR_vector, 10*log10(OMP_MSE));

% GAMP_no_quantize_MSE = mean(GAMP_no_quantize_MSE_All,1);
% plot(SNR_vector, 10*log10(GAMP_no_quantize_MSE), 'k');

GAMP_NMSE = mean(GAMP_SE_All,1)./mean(h_power_All, 1);
figure,
plot(SNR_vector, 10*log10(GAMP_NMSE), 'r');

figure
subplot(121);
bar3(abs(H_Sparse));
subplot(122);
bar3( reshape (sqrt(sum( (reshape(xhat,[], 2)).^2, 2)), Nr, Nt)) ;


% Lasso_MSE = mean(min(Lasso_MSE_All, [], 3));
% plot(SNR_vector, 10*log10(Lasso_MSE));
% plot(-10:2:10, [
%   -11.5348
%   -13.3675
%   -15.0534
%   -16.6589
%   -18.2102
%   -19.4684
%   -20.7807
%   -21.7634
%   -21.7911
%   -19.4681
%   -17.3267])
% grid on;
% 10*log10(ML_MSE)
% 10*log10(Lasso_MSE)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% sparseAWGN:  Example of estimating a sparse vector with Gaussian noise.
%
% In this problem, x is a Bernoulli-Gaussian random vector, that we want to
% estimate from measurements of the form
%
%   y = A*x + w,
%
% where A is a random matrix and w is Gaussian noise.  This is a classical
% compressed sensing problem of estimating a sparse vector for random
% linear measurements.

clc;
clear;
% close all;
Nt = 10;
Nr = 10;
L = 256;
sparseRat = 0.05;
load Data.mat

% Add to path
addpath('../main/');

SNR_vector = -10:2:20;
for index_H = 1:1:1
    index_H
    H = H_all(1:Nr, 1:Nt, index_H);
    H_Sparse = H.*(Sparse_num(1:Nr, 1:Nt) < sparseRat);
    %     H_Sparse = H_v;
        h_cv = reshape(H_Sparse, [], 1);
    %     h = [real(h_cv); imag(h_cv)];
    %     noise_cv = reshape(Noise(1:Nr, 1:L), [], 1);
    %     noise = [real(noise_cv); imag(noise_cv)];
    
    for index_SNR=1:1:length(SNR_vector)
        SNR = SNR_vector(index_SNR);
        SNR_linear = 10^(SNR/10);
        
        %         Walsh_Hadamard_matrix = hadamard(L);
        %         Z_train_scaled = sqrt(SNR_linear/Nt)*Walsh_Hadamard_matrix(1:Nt, 1:L);
        
        Z_train = (randn(Nt, L) + 1j* randn(Nt, L));
        
        Z_train_scaled = Z_train./sqrt(sum(sum(abs(Z_train).^2)))*sqrt(L)*sqrt(SNR_linear);
        
        r = reshape( dftmtx(Nr)/sqrt(Nr) * H_Sparse * Z_train_scaled + Noise(1:Nr, 1:L), [],1);
        y_train = sign( real(r)) + 1j * sign(imag(r));
        
        % Generate input estimation class
        % First, create an estimator for a Gaussian random variable (with no
        % sparsity)
        map = false;
        BGmean = 0;
        BGvar = 1;
        inputEst0 = CAwgnEstimIn(BGmean, BGvar, map, 'autoTune', true);
        
        
        % Then, create an input estimator from inputEst0 corresponding to a random
        % variable x that is zero with probability 1-sparseRat and has the
        % distribution of x in inputEst0 with probability sparseRat.
        inputEst = SparseScaEstim(inputEst0, sparseRat, 0, 'autoTune', true);
        
        % Run the GAMP algorithm
        % Set the default options
        opt = GampOpt();
        opt.nit = 50;
        opt.tol = max(min(1e-3, 10^(-SNR/10)),1e-15);
        opt.uniformVariance=0;
        opt.pvarMin=0;
        opt.xvarMin=0;
        opt.adaptStep=true;
        opt.adaptStepBethe=true;
        opt.legacyOut=false;
        
        % Create the EstimOut class object for the probit channel (MMSE)
        hA = @(in) Afast_Narrowband(in,Z_train_scaled,Nr,1);
        hAt = @(in) Afast_Narrowband(in,Z_train_scaled,Nr,0);
        lA = FxnhandleLinTrans(Nr*L,Nt*Nr,hA,hAt); % linear transform object
        probit_var = 1;
        maxSumVal = false;
        outputEst = CProbitEstimOut(y_train, 0, probit_var, maxSumVal);
        [estFin,optFin,estHist] = gampEst(inputEst, outputEst, lA, opt);
        xhat = estFin.xhat;
        GAMP_SE_All(index_H, index_SNR) = norm(xhat-h_cv)^2;
        h_power_All(index_H, index_SNR) = norm(h_cv)^2;
        %         GAMP_MSE_All(index_H, index_SNR) = (norm(xhat-h))^2/(norm(h)^2);
    end;
end;
% figure,
% ML_MSE = mean(ML_MSE_All,1);
% plot(SNR_vector, 10*log10(ML_MSE));
% hold on;

% LMMSE_MSE = mean(LMMSE_MSE_All)
% plot(SNR_vector, 10*log10(LMMSE_MSE), 'r');

% OMP_MSE = mean(OMP_MSE_All);
% plot(SNR_vector, 10*log10(OMP_MSE));

% GAMP_no_quantize_MSE = mean(GAMP_no_quantize_MSE_All,1);
% plot(SNR_vector, 10*log10(GAMP_no_quantize_MSE), 'k');

GAMP_NMSE = mean(GAMP_SE_All,1)./mean(h_power_All, 1);
figure,
plot(SNR_vector, 10*log10(GAMP_NMSE), 'r');

figure
subplot(121);
bar3(abs(H_Sparse));
subplot(122);
bar3( reshape ( abs(xhat) , Nr, Nt)) ;

% Lasso_MSE = mean(min(Lasso_MSE_All, [], 3));
% plot(SNR_vector, 10*log10(Lasso_MSE));
% plot(-10:2:10, [
%   -11.5348
%   -13.3675
%   -15.0534
%   -16.6589
%   -18.2102
%   -19.4684
%   -20.7807
%   -21.7634
%   -21.7911
%   -19.4681
%   -17.3267])
% grid on;
% 10*log10(ML_MSE)
% 10*log10(Lasso_MSE)


