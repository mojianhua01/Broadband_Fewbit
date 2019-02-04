clc;
clear all;

for i=1:1:10
    a = randn(1);
    normcdf( a )
    1/2 * erfc( - a /sqrt(2))
end;