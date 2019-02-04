function [pnSequence,chipIndex,N]=genPNSequence(polynomial,initialCondition)
%function to generate a maximal length PN sequence (m-sequence)
%Argument : Examples for polynomial input is as follows
%Polynomial [X^3+X^1+1] is given as [3 1 0]
%Polynomial [X^6+X^1+1] is given as [6 1 0]
%Generates the PN sequence,chip Index and the length of the Sequence N
%polynomial=[6 1 0]; %example: AX^3+BX^1+1 is represented as [X3 X2 X1 X0]->[A B 0]
N=2^polynomial(1)-1;
%convert polynomial representation to BnX^n+Bn-1x^n-1+...B1X^1+B0X^0 form
temp=zeros(1,max(polynomial)+1);
for i=1:length(polynomial)
    temp(length(temp)-polynomial(i))=1;
end
polynomial=temp;

if nargin<2
    initialCondition=[1 zeros(1,length(polynomial)-2)];
end
%initialize the initial state of the LFSR 
x=initialCondition;
y=zeros(1,N);

for i=1:N
    y(i)=x(end); %output;
    xi=mod(sum(x.*polynomial(2:end)),2);
    x(2:end)=x(1:end-1);
    x(1)=xi;
end
pnSequence = y;
chipIndex = 1:N;

for i=1:N
    y(i)=x(end); %output;
    xi=mod(sum(x.*polynomial(2:end)),2);
    x(2:end)=x(1:end-1);
    x(1)=xi;
end
pnSequence = y;
chipIndex = 1:N;