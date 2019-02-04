function [crossCorrelation,lags]=crossCorr(sequence1,sequence2)

%Covert to polar form
a1 = 1 - 2*sequence1;
a2 = 1 - 2*sequence2;
a1_present = a1;
a2_present = a2;
L = length(a1_present);

%Compute the cross correlation
a_past = a2_present'; %delay=0
a_present=a1_present;
for k = 1:L
    C(k) = (a_present*a_past)/L;
    a_past_out = a_past(end);
    a_past(2:end) = a_past(1:end-1);
    a_past(1) = a_past_out;
end

%Computing the simulated Cross correlation using Conjugate Symmmetry property
C_sim = [conj(fliplr(C(2:end))) C];
%delay vector
l = -(L-1):(L-1);

crossCorrelation=C_sim;
lags=l;