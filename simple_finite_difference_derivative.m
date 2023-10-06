function [deriv] = simple_finite_difference_derivative(vals, h)
% This function gives a simple approximation for the derivative
% It takes in an array of length n and a grid spacing h, and outputs an
% array of length n-1 with a finite difference approximation for the
% derivative

n = length(vals);

deriv = zeros(n-1, 1);

for i = 1:n-1
    deriv(i) = (vals(i+1) - vals(i))/h;
end

% deriv = [0; deriv];

end