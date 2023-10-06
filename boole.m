function [area] = boole(vals, h)
% This method approximates area under a curve using Boole's Rule (Deluxe
% Simpson's, exact up to cubic (?) polynomials)

% It takes in an array of function values and a mesh size, then returns an
% approximation for the integral.

area = 0;

for i = 1:4:length(vals)-4
    area = area + 2/45 * h * (7*vals(i) + 32 * vals(i+1) + 12 * vals(i+2) + 32 * vals(i + 3) + 7*vals(i+4));
end

end