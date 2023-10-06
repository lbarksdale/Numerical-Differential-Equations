function [area] = left_side_sum(vals, h)
% This function is a quadrature method using left rectangle approximation
% It takes in an array of function values and approximates the integral
% over the entire domain using left side rectangles, with a mesh sizing of h

% Note that the last function value in the domain should NOT be included in
% the input function.

area = 0;

for i = 1:length(vals)
    area = area + vals(i) * h;
end

end