%Main Program - Execute this
clc; clear all;
n=6;
N=2^n-1;
if rem(n,2) || (mod(n,4)==2)
    disp('Condition 1 Satisfied');
end
q=5;
[d,dDecimated]=genDecimatedPNSequence([6 1 0],q);

figure;
subplot(2,1,1);
plot(d,'-sr','markersize',4);
title('Preferred Pair of m-sequences');
ylabel('PN Sequence - 1');
xlabel('Chip Index (k)');
ylim([-0.2 1.2]);
grid on;

subplot(2,1,2);
plot(dDecimated,'-*k','markersize',4);
title('Preferred Pair of m-sequences');
ylabel('PN Sequence - 2');
xlabel('Chip Index (k)');
ylim([-0.2 1.2]);
grid on;

% simulated Cross Correlation
[C_sim, lags] =crossCorr (d, dDecimated);
figure;
plot(lags,C_sim,'-*r','markersize',4);
title('Cross Correlation of Preferred Pair of m-sequences');
ylabel('C_sim');
xlabel('Lag (l)');
grid on;

%Theoretical Cross Correlation
if rem(n,2)
    tn=1+2^(0.5*(n+1)); %For n=odd
else
    tn=1+2^(0.5*(n+2)); %For n=even
end
theoreticalCorr=[-1/N*tn,-1/N,1/N*(tn-2)]
