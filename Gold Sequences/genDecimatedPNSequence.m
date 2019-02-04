function [pnSequence,decimatedVersion]=genDecimatedPNSequence(polynomial,decimationFactor) 
%function to generate a decimated maximal length PN sequence (m-sequence) 
q=decimationFactor; 
pnSequence=genPNSequence(polynomial); 
dRep=repmat(pnSequence,1,q); 
decimatedVersion=dRep(1,q:q:end); 