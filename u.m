function [solution] = u(t, x, k)
% A brief function that returns the boundary data and exact solution for
% problem 6.3.10 of Strikwerda (Assignment 8 of Math 6318)

% The true solution is based on an infinite sum; this approximation cuts
% off the sum after k terms

solution = 1/2;
k = k-1;
for i = 0:k
    solution = solution + 2*(-1)^i * (cos(pi * (2*i+1) * x) / (pi * (2*i+1))) * exp(- (pi)^2 * (2*i + 1)^2 * t);
end

end