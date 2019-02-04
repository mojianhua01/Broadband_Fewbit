function frame = func_ZadoffChuSeq(Nt, N)

% if isprime(N) == 0
%     error('N should be a prime number.')
% end;

frame = zeros(Nt, N);
if isprime(N) == 1 % prime length Zadoff-Chu sequences have very good cross-correlation property
    m = (0:N-1);
    prime_num = primes(N);
    
    for i=1:1:Nt
        R = prime_num(i);
        frame(i,:) = exp( -1i * pi * R * m.*(m+1) / N );
    end;
else
    prime_num = primes(N);
    N_max = prime_num(end);
    m = (0:N_max-1);
    for i=1:1:Nt
        R = prime_num(i);
        frame(i,1:N_max) = exp( -1i * pi * R * m.*(m+1) / N_max );
        frame(i, N_max+1:1:N) = frame(i, 1:1:N-N_max);
    end;
end;