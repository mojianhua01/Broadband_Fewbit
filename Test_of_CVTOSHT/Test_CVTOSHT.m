%%function [ output_args ] = func_Broadband_Few_Bit_UPA_CVTOSHT( input_args )
%FUNC_BROADBAND_FEW_BIT_UPA_CVTOSHT Summary of this function goes here
%   Detailed explanation goes here

M = 60;
N = 100;
K = 15;
wvar = 0.05;

A = randn(M,N)/sqrt(N);

% draw signal realization
support = randperm(N,K);
x = zeros(N,1);

x(support) = (randn(K,2)*[1;1j])/sqrt(2);

% normalize signal
if 1
    xvar0 = N/K; % desired variance of non-zero coefficients
    x = sqrt(xvar0)*x; % make squared-norm ~= N
else
    x = x/norm(x)*sqrt(N); % make squared-norm == N
end

w = randn(M,2)*[1;1j]/sqrt(2);

w = sqrt(wvar)*w; % make variance = wvar;

% compute pre-quantization outputs
z = A*x; % noiseless unquantized transform outputs
y = z+w; % received, pre-quantized signal

stepsize_vector = [1.5956 0.9957 0.586 0.3352 0.1881 0.1041 0.0569 0.0308];
beta_vector = [0.3634 0.1188 0.03744 0.01154];

bit = 1;

if bit == +inf
    r = y;
else
    if bit == 1;
        stepsize = rms([real(y); imag(y)])* stepsize_vector(bit);
        %     stepsize = max(abs([real(y); imag(y)])) * 1.01;  % if ADC has only one bit, set the stepsize correctly
    else
        stepsize = rms([real(y); imag(y)])* stepsize_vector(bit);
    end;
    
    r = sign(real(y)) .* ( min( ceil( abs(real(y)) /stepsize) , 2^bit-1 ) - 1/2 ) * stepsize  + ...
        1j* sign(imag(y)) .* ( min( ceil( abs(imag(y)) /stepsize) , 2^bit-1 ) - 1/2 ) * stepsize;
    
end;


% Yaniv cross-validation
Kxval_folds = 3;
Kxval_grid_length = 10;
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
tstart = tic;
if Kxval_grid_length>1
    
    err = nan(Kxval_folds,Kxval_grid_length); % test error counts
    m = randperm(M); % randomly ordered output indices
    for f=1:Kxval_folds,
        % set test/train indices
        i_min = 1+round(M*(f-1)/Kxval_folds);
        i_max = round(M*f/Kxval_folds);
        m_test = m(i_min:i_max); % testing indices
        m_train = m([1:i_min-1,i_max+1:M]); % training indices
        
        % evaluate performance over grid
        xhat = (A(m_train,:).'*r(m_train)); % non-sparse training estimate
        [~,indx] = sort(abs(xhat),'descend');
        for k=1:Kxval_grid_length
            % estimate for Khat sparsity hypothesis
            Khat = Kgrid_xval(k);
            xhatY = zeros(N,1);
            xhatY(indx(1:Khat)) = xhat(indx(1:Khat)); % sparse training estimate
            
            % evaluate performance on test indices
            y_test = A(m_test,:)*xhatY;
            
            r_test = sign(real(y_test)) .* ( min( ceil( abs(real(y_test)) /stepsize) , 2^bit-1 ) - 1/2 ) * stepsize  + ...
                1j* sign(imag(y_test)) .* ( min( ceil( abs(imag(y_test)) /stepsize) , 2^bit-1 ) - 1/2 ) * stepsize;
            err(f,k) = norm(r_test - r(m_test)).^2;
        end
    end
    
    % best sparsity estimate
    [~,kxval] = min(sum(err,1));
    Kxval = Kgrid_xval(kxval);
    
else % bypass cross-validation and use single value in grid
    
    Kxval = Kgrid_xval;
    
end

% final estimation
xhat = (A.'*r); % non-sparse estimate
[~,indx] = sort(abs(xhat),'descend');
xhatY = zeros(N,1);
xhatY(indx(1:Kxval)) = xhat(indx(1:Kxval)); % sparse estimate

figure, stem(abs(xhatY))
figure, stem(abs(x))

% scale based on power estimate
xhatY = xhatY*(sqrt(N*Ez2hat)/norm(xhatY));
runtimeY = toc(tstart);

% evaluate output
if knowTrue
    nmse_xY = (norm(xhatY-xTrue)/norm(xTrue))^2;
    nmse_xY_scaled = norm(xhatY/norm(xhatY)-xTrue/norm(xTrue))^2;
end
sparsityFinY = Kxval/N;


% end

