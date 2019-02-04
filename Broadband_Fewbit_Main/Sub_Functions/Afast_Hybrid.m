% out = Afast(in,B,Nr,mode)
%
% When mode==0, fast implementation of multiplication by matrix 
%     A = kron(transpose(B), dftmtx(Nr)/sqrt(Nr));  % A_train_cv is complex-valued
%     z = A *x;      % "True" transform vector
% When mode~=0, fast multiplication by its Hermitian.
%
function out = Afast_Hybrid(in, Analog_combining_pattern, B, Nr, mode)

%  Nr = size(in,1);
%  T = size(B,2); 
%  N = size(B,1);
%  Nt = size(B,3);

 if (mode==1) % perform regular multiplication

%    out = kron(transpose(B), dftmtx(Nr)/sqrt(Nr)) * in;
    temp = dftmtx(Nr)/sqrt(Nr) * reshape(in, Nr, []) * B;
    out = reshape( temp.*Analog_combining_pattern, [], 1);

 else % perform Hermitian transpose

%    out = kron(transpose(B), dftmtx(Nr)/sqrt(Nr))' * in;
    temp = dftmtx(Nr)'/sqrt(Nr) * ( reshape(in, Nr, []).*Analog_combining_pattern) * B';
    out = reshape( temp, [], 1);
 end

