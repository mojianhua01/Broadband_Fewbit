% evaluates channel estimate in terms of mutual info after various types of precoding

function [MI, Rate] = eval_chanest(Ghat,Gtrue, Params)

% extract dimensions
Nt = Params.Nt;
Nr = Params.Nr;
Nd = Params.Nd;
Bt = Params.Bt;
Br = Params.Br;
N_sc = Params.N_sc;
SNR_linear = Params.SNR_linear;
Nmin = min(Nt,Nr);
beta = Params.beta;
Nc = Params.Nc;
Np = Params.Np;

%  % convert to frequency domain and perform noise-whitening
%  FHhat = bsxfun(@times,reshape(sqrt(Nr/wvar)./winr',[1,Nr,1]),...
%                         ifft(fft(Ghat,N,3),Nr,2)); % Nt x Nr x N
%  FHtrue = bsxfun(@times,reshape(sqrt(Nr/wvar)./winr',[1,Nr,1]),...
%                         ifft(fft(Gtrue,N,3),Nr,2)); % Nt x Nr x N

% From angular domain to spatial domain
% From spatial domain to frequency domain
% FD: frequency domain
for index_delay = 1:1:size(Gtrue,3)
    Htrue(:, :, index_delay) = Br * ...
        Gtrue(:,:, index_delay) * Bt';
end;
FHtrue = fft(Htrue, N_sc, 3);

for index_delay = 1:1:size(Gtrue,3)
    Hhat(:, :, index_delay) = Br * ...
        Ghat(:,:, index_delay) * Bt';
end;
FHhat = fft(Hhat, N_sc, 3);  

% SVD at each subcarrier
Sighat = zeros(Nmin,N_sc);
Uhat = cell(1,N_sc);
Vhat = cell(1,N_sc);
Sigtrue = cell(1,N_sc);
Utrue = cell(1,N_sc);
Vtrue = cell(1,N_sc);
for n=1:N_sc
    [Uhat{n},tmp,Vhat{n}] = svd(FHhat(:,:,n),'econ');
    Sighat(:,n) = diag(tmp);
    [Utrue{n},Sigtrue{n},Vtrue{n}] = svd(FHtrue(:,:,n),'econ');
end

% waterfilling: argmax_p sum(log(1+p.*snrhat)) s.t. sum(p)=Pmax & p>=0
snrhat = (Sighat(:).^2).';   

Pmax = SNR_linear * N_sc;  % the total signal power is normalized by N_sc if the noise power is kept fixed
% the power is multipled by N_sc because there are N_sc symbols in the time-domain
% the noise power should be normalized by 1/N_sc if the total signal power is kept fixed

PowerAllo = (Pmax + sum(1./snrhat))/length(snrhat) - 1./snrhat;
while(~isempty(find(PowerAllo < 0, 1 )))
    IndexN = (PowerAllo <= 0) ;
    IndexP = find(PowerAllo > 0);
    MP = length(IndexP);
    PowerAllo(IndexN) = 0;
    ST = snrhat(IndexP);
    PowerAlloT = (Pmax + sum(1./ST))/MP - 1./ST;
    PowerAllo(IndexP) = PowerAlloT;
end

p = reshape(PowerAllo, Nmin, N_sc);

%% joint decoding
% bpcu1 = zeros(1,N_sc);
% for n=1:N_sc
%     Rx = Vhat{n} * diag(p(:,n)) * Vhat{n}';
%     Rn = (1-beta) * eye(Nr) + beta *(1-beta)* diag(diag(FHtrue(:,:,n) *Rx* FHtrue(:,:,n)'));
%     bpcu1(n) = real( 1/N_sc * log2( det(eye(Nr) + (1-beta)^2* Rn^(-1)*FHtrue(:,:,n)*Rx* FHtrue(:,:,n)')) );
% end;
% rate1 = sum(real(bpcu1));

% Rx = zeros(Nt, Nt);
% for n=1:N_sc
%     Rx = Rx + 1/N_sc * Vhat{n} * diag(p(:,n)) * Vhat{n}';
% end;
% %% each stream is decoded independently
% for n=1:N_sc
% %     Rx = Vhat{n} * diag(p(:,n)) * Vhat{n}';
%     Rn = (1-beta) * eye(Nr) + beta *(1-beta)* diag(diag(FHtrue(:,:,n) *Rx* FHtrue(:,:,n)'));
%     Rnn = Uhat{n}' * Rn * Uhat{n};
%     Gamn = Uhat{n}'*FHtrue(:,:,n)*Vhat{n}*sqrt(diag((p(:,n))));
%     for m = 1:1:Nmin
%         bpcu(n,m) = real( 1/N_sc*log2( 1 + (1-beta)^2* abs(Gamn(m,m))^2 *...
%             ( Rnn(m,m) + sum(abs(Gamn(m,:)).^2) - abs(Gamn(m,m))^2 )^(-1) )); % rate of the m-th stream of the n-subcaarrier
%     end;
% end;
% 
% rate = sum(sum(real(bpcu)));

Rx = cell(1, N_sc);
for n=1:N_sc
    Rx{n} = Vhat{n} * diag(p(:,n)) * Vhat{n}';
end;
%% each stream is decoded independently
for n=1:N_sc
%     Rx = Vhat{n} * diag(p(:,n)) * Vhat{n}';
    Rn = (1-beta) * eye(Nr) + beta *(1-beta)* diag(diag(FHtrue(:,:,n) *Rx{n}* FHtrue(:,:,n)'));
    Rnn = Uhat{n}' * Rn * Uhat{n};
    Gamn = Uhat{n}'*FHtrue(:,:,n)*Vhat{n}*sqrt(diag((p(:,n))));
    for m = 1:1:Nmin
        bpcu(n,m) = real( 1/N_sc*log2( 1 + (1-beta)^2* abs(Gamn(m,m))^2/...
            ( Rnn(m,m) + (1-beta)^2 * (sum(abs(Gamn(m,:)).^2) - abs(Gamn(m,m))^2 )) )); % rate of the m-th stream of the n-subcaarrier
    end;
end;

MI = sum(sum(real(bpcu)));
Rate = (Nc - Np)/Nc * MI;

return;
