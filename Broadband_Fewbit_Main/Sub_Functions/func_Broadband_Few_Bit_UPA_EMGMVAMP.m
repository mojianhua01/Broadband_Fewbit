function [ ghat, time_EM_GM_VAMP] = func_Broadband_Few_Bit_UPA_EMGMVAMP( T_cyclic, r_bit, Params)
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

Nc = Params.Nc;
Np = Params.Np;

bit = Params.bit;
stepsize = Params.stepsize;
SNR_linear = Params.SNR_linear;
dither_mean = Params.dither_mean;
wvar = Params.wvar;

Gtrue = Params.Gtrue;

% Bt = Params.Bt;
% V_cyclic_0 = kron(eye(Nd), Bt') * T_cyclic;
% fast implementation of V_cyclic = kron(eye(Nd), Bt') * T_cyclic
% T_temp = reshape(T_cyclic, Nt, Nd, Np);
% V_temp = Bhfast(T_temp, Nt);
% V_cyclic = reshape(V_temp, Nd*Nt, Np);
V_cyclic = Params.V_cyclic;

% % EG-BG-GAMP setup
sparseRat = 0.05;    % Bernoulli-Gaussian active probability
BGmean = 0;        	% Bernoulli-Gaussian active mean
BGvar = 1/sparseRat/Nd;         % Bernoulli-Gaussian active variance

% Decide on MAP or MMSE GAMP. Only support MMSE GAMP in our codes!
map = false;

% Auto tune the mean, variance and sparseRat.
% This is extremely important for the algorithm to work!!!

%% EG-GM-VAMP setup
if 1
    Lgm = 3;            % number of GM components
    omega0 = ones(Nt*Nr*Nd, 1, Lgm)/Lgm;
    %theta0 = bsxfun(@times, ones(Nt*Nr*Nd, 1, Lgm),...
    %                reshape(linspace(-1,1,Lgm)/sqrt(1/3),[1,1,Lgm]) );
    %phi0 = ones(Nt*Nr*Nd, 1, Lgm)/Lgm;
    theta0 = BGmean*ones(Nt*Nr*Nd, 1, Lgm);
    vars = [1:Lgm]/sqrt(Lgm);
    phi0 = bsxfun(@times, ones(Nt*Nr*Nd, 1, Lgm),...
        reshape(BGvar*vars/mean(vars),[1,1,Lgm]) );
    inputEst0 = CGMEstimIn(omega0, theta0, phi0, 'autoTune', true, 'thetaTune', false);
else % BG case for debugging
    inputEst0 = CAwgnEstimIn(BGmean, BGvar, map, 'autoTune', true, 'mean0Tune', false);
end
inputEst = SparseScaEstim(inputEst0, sparseRat, 0, 'autoTune', true);


if bit == +inf
    outputEst = CAwgnEstimOut(r_bit, wvar, map);
else % uniform quantization
    % Create an output estimation class corresponding to the Gaussian noise.
    % Note that the observation vector is passed to the class constructor, not
    % the gampEst function.
    % outputEst = CQuantizeEstimOut(r_bit, bit, stepsize, 0, wvar, map);
    outputEst = CQuantizeEstimOut(r_bit, bit, stepsize, dither_mean, wvar, map);
end;

%% Set the default options
if bit == +inf
    maxIt = 50;
    tol = max(min(1e-3, 1/SNR_linear),1e-15);
else
    maxIt = 25; %25
    tol = 10^(-2.5); % 10^(-2.5)
end
damp = 0.9;

%% setup VAMP
vampOpt = VampGlmOpt;
vampOpt.nitMax = maxIt;
vampOpt.tol = tol;
vampOpt.damp = damp;
vampOpt.verbose = false;
vampOpt.debug = false;
vampOpt.fxnErr1 = @(x1,z1) 10*log10( sum(abs(x1-Gtrue(:)).^2,1)./sum(abs(Gtrue(:)).^2,1) );
vampOpt.fxnErr2 = @(x1,z1) 10*log10( sum(abs(...
    bsxfun(@times, x1, sum(conj(x1).*Gtrue(:),1)./sum(abs(x1).^2,1)) - Gtrue(:)...
    ).^2,1)./sum(abs(Gtrue(:)).^2,1) );


%% fast implementation of matrix-vector multiplication

%% UPA
hA = @(in) Afast_Digital_UPA(in, T_cyclic, V_cyclic, Params, 1);
hAt = @(in) Afast_Digital_UPA(in, T_cyclic, V_cyclic, Params, 0);

%lA = FxnhandleLinTrans(Nr*Np, Nt*Nr*Nd, hA, hAt); % linear transform object
% A has size Nr*Np  X  Nt*Nr*Nd

% There's something important I forgot to say about VAMP.  The algorithm uses an SVD to avoid the need to do a matrix inverse at each iteration.  So for VAMP to enjoy the computational benefits of a fast linear operator, we need that the SVD operations are themselves fast.
% To explain what I mean, let's suppose that A = U*S*V' defines the linear operator and its SVD.  Given matrices {U,S,V}, the VAMP algorithm can be written in various ways that involve vector-multiplications with the matrices U and V (and in some cases also A).
% But sometimes we want to represent A as a function handle rather than as an explicit matrix, one that calls a fast algorithm like an FFT.  This is precisely what happened in our mmW work.  So, to make VAMP fast, we need function handles for U and V that call fast algorithms.
% Fortunately this occurs in the case of ZC sequences.  There, A is a tall matrix with orthogonal columns.  Thus the SVD takes the form A = U*S*V' with V=I, S=[c*I;0], some scalar c, and some unitary U.  But we would like to avoid explicitly deriving/using the large square U matrix.  Fortunately it turns out that VAMP is coded such that, when A is tall, it only needs function handles to A and V.  Likewise, when A is wide, it only needs function handles to A and U.
% The upshot is that, when you call VAMP for the ZC case, you only need to supply it with a function handle to A (which you already have from the GAMP case) and the trivial function handle "V = @(x) x;".  You will also need to supply it with the eigenvalues of A'*A, which you know in closed form: they are the square of the singular values that we specify after (42) in our paper.
% For the other training sequences, I recommend explicitly computing the eigenvalue decomposition "[V,D]=eig(A'*A)" offline and saving the V,D quantities in a file.  (This may take a while depending on the dimensions!)  For numerical reasons, you may actually want to do "AhA=A'*A; [V,D]=eig(0.5*(AhA+conj(AhA)))".

vampOpt.Ah = hAt;
vampOpt.gam1xinit = min(1, 1e1/(SNR_linear)); %1e-2
vampOpt.gam1zinit = min(1, 1e1/(SNR_linear)); %1e-2
vampOpt.N = Nt*Nr*Nd;
if strcmp(Params.Sequences, 'Shifted_ZC')
    
    vampOpt.d = (SNR_linear * Np/Nt) * ones(Nt*Nr*Nd, 1);
    vampOpt.V = @(x) x;
    vampOpt.Vh = @(x) x;
    
    %% stop function: stop when the NMSE is smaller than that of BPDN (i.e., SPGL1) algorithm
%     NMSE_BPDN_dB_vector=[-11.67099606	-12.57338142	-12.91700443	-13.77483569	...
%         -14.29506195	-14.66714104	-14.95196809	-15.34087005	-15.5699147]; % 0dB 4-bit shifted-ZC
%     Np_vector = 512*(2:1:10);
%     NMSE_BPDN_dB = NMSE_BPDN_dB_vector( Np_vector==Np );
%     vampOpt.fxnStop = @(i,err1,err2,...
%         r1old,r1,gam1x,x1,eta1x,...
%         p1old,p1,gam1z,z1,eta1z,...
%         r2old,r2,gam2x,x2,eta2x,...
%         p2old,p2,gam2z,z2,eta2z) (err1<NMSE_BPDN_dB); % err1<-10
else
    
    %     vampOpt.d = Params.D_of_A.^2;
    %     vampOpt.V = @(x) Params.V_of_A * x;
    %     vampOpt.Vh = @(x) Params.V_of_A' * x;
    
    V_of_VcVt = Params.V_of_VcVt;
    D_of_VcVt = Params.D_of_VcVt;
    vampOpt.d = kron(D_of_VcVt, ones(Nr, 1));
    vampOpt.V = @(x) reshape ( reshape(x, Nr, Nd*Nt) * transpose(V_of_VcVt),  [], 1); % @(x) kron(V_of_V_cyclic, eye(Nt)) * x = eye(Nt) * X * transpose(V_of_V_cyclic)
    vampOpt.Vh = @(x) reshape ( reshape(x, Nr, Nd*Nt) * conj(V_of_VcVt) ,  [], 1);
    
%     vampOpt.fxnStop = @(i,err1,err2,...
%         r1old,r1,gam1x,x1,eta1x,...
%         p1old,p1,gam1z,z1,eta1z,...
%         r2old,r2,gam2x,x2,eta2x,...
%         p2old,p2,gam2z,z2,eta2z) (err1<-80);  %err1<-10
end;




%% Run the VAMP algorithm
mode = 1;
switch mode
    case 1
        tic;
        Reps = 1;
        for i=1:1:Reps
            ghat = VampGlmEst(inputEst, outputEst, hA, vampOpt);
        end;
        time_EM_GM_VAMP = toc/Reps
    case 2  % to plot the figure of performance-time trade-off
        for iteration_num = 1:1:13
            vampOpt.nitMax = iteration_num+1;
            vampOpt.tol = 0;
            tic;
            Reps = 10;
            for i=1:1:Reps
                ghat = VampGlmEst(inputEst, outputEst, hA, vampOpt);
            end;
            time_EM_GM_VAMP(iteration_num) = toc/Reps
        end;
    case 3 % debug mode
        tic;
        vampOpt.nitMax = 51;
        vampOpt.tol = 0;
        [ghat,estFin] = VampGlmEst(inputEst, outputEst, hA, vampOpt);
        time_EM_GM_VAMP = toc
        
        figure,
        plot(estFin.err1)
        estFin.err1'
        xlabel('Iteration')
        ylabel('Normalized MSE')
        title('EM-GM-VAMP')
        grid on;
end
end

