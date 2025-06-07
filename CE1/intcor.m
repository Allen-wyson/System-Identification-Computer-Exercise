function [R,h] = intcor(u,y)
%INTCOR Summary of this function goes here
%   Detailed explanation goes here

M = length(u); % M est la période du signal
h = -floor(M/2):floor(M/2);
R = zeros(1, length(h)); 


for i = 1:length(h)
    sum_value = 0;

    for k = 1:M
        y_shifted = y(mod(k - h(i) - 1, M) + 1);
        sum_value = sum_value + u(k) * y_shifted;
    end
    
    R(i) = sum_value / M; 
end
