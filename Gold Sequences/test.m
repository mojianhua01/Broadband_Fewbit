clc;
close all;
% dbstop if error
P = 509;
N = 509;
seq = lteZadoffChuSeq(5,P);
seq = [seq ; seq(1:N-P, 1) ];
figure,
plot(abs(fcxcorr(seq, seq)./length(seq)));

seq2 = lteZadoffChuSeq(7,P);
seq2 = [seq2 ; seq2(1:N-P, 1) ];
figure,
plot(abs(fcxcorr(seq2, seq2)./length(seq2)));

out = fcxcorr(seq,seq2);

figure, plot(abs(out));
title('cross correlation')

primes(P);
isprime(P);

seq3 = randn(N,1);
figure,
plot(abs(xcorr(seq3)./length(seq3)));
seq4 = randn(N,1);
out = fcxcorr(seq3,seq4);

figure, plot(abs(out));
title('random sequences')


goldseq1 = comm.GoldSequence('FirstPolynomial','x^5+x^2+1',...
    'SecondPolynomial','x^5+x^4+x^3+x^2+1',...
    'FirstInitialConditions',[0 0 0 0 1],...
    'SecondInitialConditions',[0 0 0 0 1],...
    'Index',1,'SamplesPerFrame',31);
x1 = step(goldseq1)

goldseq2 = comm.GoldSequence('FirstPolynomial','x^5+x^2+1',...
    'SecondPolynomial','x^5+x^4+x^3+x^2+1',...
    'FirstInitialConditions',[0 0 0 0 1],...
    'SecondInitialConditions',[0 0 0 0 1],...
    'Index',2,'SamplesPerFrame',31);
x2 = step(goldseq2)

out = fcxcorr(2*x1-1,2*x2-1);

figure, plot((out));


% clear all
% 
% n=input('please input the number:')                     %??????
% 
% m=2;                                                                %??????2??
% display([num2str(n),'='])
% 
% while(1)
%       if(~mod(n,m))                                              %????????
%               k=m;
%               if(n==k)                                               %????????
%                     display([num2str(n)])
%                     break;                                            %????
%               else
%                     n=n/k;                                           %?n??????????
%                     m=1;                                             %???????2??
%                     display([num2str(k),'*'])                  %??????????
%               end
%       end
%      
% 
%       m=m+1;                                  %????m??????????4?6?????????2?3????????
% end