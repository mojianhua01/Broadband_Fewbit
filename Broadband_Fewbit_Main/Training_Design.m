% clear;
% close all;
% 
% % Model parameters
% Nt = 4;
% Nr = 4;
% Nd = 10;
% srrc = false; % use SRRC instead of RC time-domain pulse?
% alf_rc = 0.5; % raised-cosine parameter [dflt=0.5]
% Npre = min(2,Nd/2); % channel taps before first arrival [dflt=5]
% Npst = min(2,Nd/2); % channel taps after last arrival [dflt=5]
% Num_cluster = 2;   % number of path
% plot_datacube = 1;
% alpha = 10;
% sparseRat = alpha * Num_cluster/Nd/Nt/Nr;
% 
% % Pilot length and SNR
% N = 2^7;
% SNR_vector = 10;
% 
% 
% BGmean = 0;        	% Bernoulli-Gaussian active mean
% BGvar = 1;         	% Bernoulli-Gaussian active variance
% 
% % create aperture/delay channel
% H = zeros(Nr,Nt,Nd); % MIMO matrix for each delay
% delay = zeros(Num_cluster,1);
% gain = zeros(Num_cluster,1);
% theta_t = zeros(Num_cluster,1);
% theta_r = zeros(Num_cluster,1);
% n = [0:Nd-1]+eps;               % baud-normalized sampling times
% for l=1:Num_cluster % for each path...
%     delay(l) = Npre+rand(1)*(Nd-Npre-Npst); % delay (between Npre and Nd-Npst)
%     if srrc
%         pulse = ( (1-alf_rc)*sinc((n-delay(l))*(1-alf_rc)) ...
%             + cos(pi*(n-delay(l))*(1+alf_rc))*4*alf_rc/pi ...
%             )./(1-(4*alf_rc*(n-delay(l))).^2);
%     else % rc pulse
%         pulse = cos(pi*alf_rc*(n-delay(l)))./(1-(2*alf_rc*(n-delay(l))).^2) ...
%             .*sinc(n-delay(l)); % raised-cosine impulse response
%     end
%     pulse = pulse/norm(pulse);
%     %plot(pulse,'.-'); title('pulse'); pause;
%     
%     theta_t(l) = rand(1); % transmit angle (between 0 and 1)
%     theta_r(l) = rand(1); % receive angle (between 0 and 1)
%     gain(l) = randn(1,2)*[1;1i]/sqrt(2*Num_cluster); % complex gain
%     for d=1:Nd % for each delay ...
%         H(:,:,d) =  H(:,:,d) + pulse(d)*gain(l)*sqrt(1/Nt/Nr) ...
%             *exp(1i*2*pi*theta_r(l)*(0:Nr-1))' ...
%             *exp(1i*2*pi*theta_t(l)*(0:Nt-1));
%         Gtrue(:,:,d) = dftmtx(Nr)'/sqrt(Nr) * H(:,:,d) * dftmtx(Nt)/sqrt(Nt);
%     end
% end %l
% 
% gtrue = reshape(Gtrue, [], 1);
% 
% noise = (randn(Nr, Nt, N) + 1i * randn(Nr, Nt, N))/sqrt(2);
% 
% 
% Walsh_Hadamard_matrix = hadamard(N);
% T_train = Walsh_Hadamard_matrix(randsample(1:N,Nt), 1:N);
% 
% % T_train = (randn(Nt, N) + 1j * randn(Nt, N));
% T_train = sign(randn(Nt,N))/sqrt(2)+1i*sign(randn(Nt,N))/sqrt(2);
% 
% % T_train = dftmtx(Nt)/sqrt(Nt)'*T_train;
% 
% for index_SNR=1:1:length(SNR_vector)
%     SNR = SNR_vector(index_SNR);
%     SNR_linear = 10^(SNR/10);
%     
%     Z_train_scaled = T_train./norm(T_train, 'fro') *sqrt(N)*sqrt(SNR_linear);
%     % scale Z_train such that sum((Z_train_scaled(:,ii)).^2) = SNR_linear for
%     % ii=1,2,...
%     
%     % create circular-delay matrices: J{d} rotates down by (d-1) places
%     J = cell(1,Nd-1);
%     J{1} = speye(N);
%     for d=1:Nd-1
%         J{d+1} = [spalloc(N-d,d,0), speye(N-d); speye(d), spalloc(d,N-d,0)];
%     end
%     
%     V = zeros(Nt*Nd, N);
%     V(1: Nt,:) = Z_train_scaled;
%     for d=2:1:Nd
%         V( (d-1)*Nt+1: d*Nt, :) = Z_train_scaled * J{d};
%     end;
%     
%     % A = kron(transpose(S), dftmtx(Nr)/sqrt(Nr));
%     % z = A *gtrue;      % "True" transform vector
%     
%     A = kron( V.', dftmtx(Nr)/sqrt(Nr) );
% end;
% 
% C = A'*A;
% figure(1), surf(abs(C))
% 
% figure(2), plot(eig(C))
% 
% D = Nt/SNR_linear/N* V*V';
% figure(3), bar3(abs(D))
% temp = D - diag(diag(D));
% var(temp(:))
% 
% figure(4), plot(eig(D))
% 
% max(max(abs(D - diag(diag(D)))/mean(diag(D))))

clear;
close all;

Np = 1024;
Nt = 32;
if floor(log2(Np)) == log2(Np)
    [a, b] = generate_golay(log2(Np)-1);
else
    keyboard;
end;
s = [a, b];
T0 = zeros(Nt,Np);
for i=0:Nt-1
    T0(i+1,:) = circshift(s,[0 i*32]);
    %                     T0(i+1,:) = [a(Np/2-i*Nd+1: Np/2), a(1:Np/2-i*Nd), ...
    %                         b(Np/2-i*Nd+1: Np/2), b(1:Np/2-i*Nd)]; %% for
    %                         this kind of T0, the singular values are same
end
[U S V] = svd(T0, 'econ');
figure, plot(diag(S))

T0_product = T0*T0';
T0_product_nondiag = T0_product - diag(diag(T0_product));
max(T0_product_nondiag(:))
min(T0_product_nondiag(:))
% figure, bar(T0*T0');