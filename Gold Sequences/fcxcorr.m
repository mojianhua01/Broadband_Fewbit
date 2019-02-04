function [ xc ] = fcxcorr(u1,u2)
%[ xc ] = fcxcorr(u1,u2)
%Uses fft to calculate the circular cross correlation of two periodic
%signal vectors.This is equivalent to xc(k)=sum(u1.*circshift(u2,k)), but
%much faster, especially for large vectors. There is no input checking; 
%vectors must be equally sized.
%The result is not normalized.  You can get the normalized result using:
% xc_n=fcxcorr(u1,u2)/(norm(u1)*norm(u2));

%copyright Travis Wiens 2009

xc=ifft(fft(u1).*conj(fft(u2)));