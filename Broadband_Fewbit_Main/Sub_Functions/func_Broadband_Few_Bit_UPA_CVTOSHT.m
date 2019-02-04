function [ xhatY ] = func_Broadband_Few_Bit_UPA_CVTOSHT( T_cyclic, r_bit, Params)
% Need revision!!!
% need to incorporate the fast implementation!
% The matrix A is so large for simulation

%FUNC_BROADBAND_FEW_BIT_UPA_CVTOSHT Summary of this function goes here
%   Detailed explanation goes here

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
Nb = Params.Nb;

Bt = Params.Bt;
Br = Params.Br;

bit = Params.bit;
stepsize = Params.stepsize;
SNR_linear = Params.SNR_linear;

% V_cyclic_0 = kron(eye(Nd), Bt') * T_cyclic;
% fast implementation of V_cyclic = kron(eye(Nd), Bt') * T_cyclic
T_temp = reshape(T_cyclic, Nt, Nd, Nb);
V_temp = Bhfast(T_temp, Nt);
V_cyclic = reshape(V_temp, Nd*Nt, Nb); 

A = kron(V_cyclic.', Br);

if bit == +inf
    r = r_bit;
else
    r = (- 2^( bit -1 ) + real(r_bit) + 1/2 ) * stepsize + ...
        1j * (- 2^( bit -1 ) + imag(r_bit) + 1/2 ) * stepsize;
end;

% Yaniv cross-validation
% [M, N] = size(A);
M = Nb*Nr;
N = Nt*Nr*Nd;

Kxval_folds = 2;
Kxval_grid_length = 20;
Kxval_grid = NaN;

% create tuning grid if not supplied by user
if isnan(Kxval_grid)
    if 0 % linearly spaced
        Kxval_grid_min = ceil(N/Kxval_grid_length);
        Kxval_grid_max = N;
        Kgrid_xval = unique(round(linspace(Kxval_grid_min,...
            Kxval_grid_max,Kxval_grid_length)));
    else % logarithmically spaced
        Kxval_grid_min = 2;
        Kxval_grid_max = N;
        Kgrid_xval = unique(round(logspace(log10(Kxval_grid_min),...
            log10(Kxval_grid_max),Kxval_grid_length)));
    end
end
Kxval_grid_length = length(Kgrid_xval);

% cross-validation tuning of K
% if Kxval_grid_length>1
%     
%     err = nan(Kxval_folds,Kxval_grid_length); % test error counts
%     m = randperm(M); % randomly ordered output indices
%     for f=1:Kxval_folds,
%         % set test/train indices
%         i_min = 1+round(M*(f-1)/Kxval_folds);
%         i_max = round(M*f/Kxval_folds);
%         m_test = m(i_min:i_max); % testing indices
%         m_train = m([1:i_min-1,i_max+1:M]); % training indices
%         
%         % evaluate performance over grid
%         xhat = (SNR_linear * Nb/Nt)^(-1) * (A(m_train,:)' * r(m_train)); % non-sparse training estimate
% %         xhat = (SNR_linear * Nb/Nt)^(-1) * (A(m_train,:)' * r(m_train)); % non-sparse training estimate
%         [~,indx] = sort(abs(xhat),'descend');
%         for k=1:Kxval_grid_length
%             % estimate for Khat sparsity hypothesis
%             Khat = Kgrid_xval(k);
%             xhatY = zeros(N,1);
%             xhatY(indx(1:Khat)) = xhat(indx(1:Khat)); % sparse training estimate
%             
%             % evaluate performance on test indices
%             y_test = A(m_test,:) * xhatY;
%             
% %             if bit == +inf
% %                 r_test = y_test;
% %             else
% %                 r_test = sign(real(y_test)) .* ( min( ceil( abs(real(y_test)) /stepsize) , 2^bit-1 ) - 1/2 ) * stepsize  + ...
% %                     1j* sign(imag(y_test)) .* ( min( ceil( abs(imag(y_test)) /stepsize) , 2^bit-1 ) - 1/2 ) * stepsize;
% %             end;
%             
%             err(f,k) = norm(y_test - r(m_test))^2;
%         end
%     end
%     
%     % best sparsity estimate
%     [~,kxval] = min(sum(err,1));
%     Kxval = Kgrid_xval(kxval);
%     
% else % bypass cross-validation and use single value in grid
%     
%     Kxval = Kgrid_xval;
%     
% end

Lloyd_stepsize_vector = [1.5956 0.9957 0.586 0.3352 0.1881 0.1041 0.0569 0.0308];
beta_vector = [0.3634 0.1188 0.03744 0.01154];

if bit == +inf;
    beta = 0;
else
    beta = beta_vector(bit);
end;

Kxval_fraction = sqrt(6)/pi *Params.estimated_norm  * ...
    sqrt(SNR_linear * Nb/ (Nt * (1-beta + beta *(1-beta) * SNR_linear) ) ) * (1-beta);
Kxval = min( round(Kxval_fraction), N);

% final estimation
xhat = (SNR_linear * Nb/Nt)^(-1) * (A'*r); % non-sparse estimate
[~,indx] = sort(abs(xhat),'descend');
xhatY = zeros(N,1);
xhatY(indx(1:Kxval)) = xhat(indx(1:Kxval)); % sparse estimate

% scale based on power estimate
xhatY = xhatY*(sqrt(Nt*Nr))/norm(xhatY);

% Kxval;
% sparsityFinY = Kxval/N;

end

