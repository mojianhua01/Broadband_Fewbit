clc;
clearvars;
y_temp = randn(10000,1);

bit_vector = [ 1 2 3 4 5 6 7 8 +inf];
Lloyd_stepsize_vector = [1.5956 0.9957 0.586 0.3352 0.1881 0.1041 0.0569 0.0308];


y_temp = randn(1e7,1);
for i = 1:1:length(bit_vector)
    bit =  bit_vector(i);
    if bit == +inf
        mu(i) = 1;
    else
        stepsize_temp =  Lloyd_stepsize_vector(bit);
        r_temp = sign( y_temp ) .* ( min( ceil( abs( y_temp )/stepsize_temp) , 2^(bit-1) ) - 1/2 ) * stepsize_temp;
        mu(i) = mean(y_temp.*r_temp)/rms(y_temp)^2;
    end;
end;

beta = 1 - mu
