function [a,b] = generate_golay(N)
% [a,b] = generate_golay(N)
%
% Generate the Golay codes a and b with length 2^N.
%
% Then write them to disk as golayA.wav and golayB.wav.


% These initial a and b values are Golay
a = [1 1];
b = [1 -1];

% Iterate to create a longer Golay sequence
while (N>1)
    olda = a;
    oldb = b;
    a = [olda oldb];
    b = [olda -oldb];

    N = N - 1;
end

% % Guess the sampling rate. It doesn't really matter.
% fs = 44100;
% 
% % Scaling by 0.9999 suppresses a warning message about clipping.
% wavwrite(a*0.9999,fs,'golayA.wav');
% wavwrite(b*0.9999,fs,'golayB.wav');