function [v] = lax_wendroff(h)
% This is a numerical solution to the one-way wave equation using the
% lax-wendroff scheme. It is designed to demonstrate the order of accuracy
% of the lax-wendroff scheme.

syms x;

% Constants
x_min = -1;
x_max = 1;
t_min = 0;
t_max = 1.2;
a = 1;
lambda = 0.8;
k = lambda * h;

% Initial data
u_0(x) = sin(2 * pi * x);

% Create array to store function values and fill lowest time level with
% initial data. Array consists of two levels, an old and a new
width = (x_max - x_min) / h;
height = (t_max - t_min) / k;
v = zeros(ceil(width)+1, 2);
for i = 1 : ceil(width)+1
    v(i,1) = u_0(x_min + h * (i-1));
end

% Old and new track which layer of the array v holds old data vs new data
old = 1;
new = 2;
% Compute next layer of time values until desired time level is reached
for j = 1 : ceil(height)+1
    % Compute all space values at given time before moving on
    % Writes new data over old (unneeded) data
    for i = 2 : ceil(width)
        v(i, new) = v(i, old) - a*lambda / 2 * (v(i+1, old) - v(i-1, old))+ a^2 * lambda^2 / 2 * (v(i+1, old) -2*v(i, old) + v(i-1, old));
    end

    % Adds periodicity
    v(ceil(width)+1, new) = v(ceil(width)+1, old) - a*lambda / 2 * (v(1, old) - v(ceil(width), old))+ a^2 * lambda^2 / 2 * (v(1, old) -2*v(ceil(width)+1, old) + v(ceil(width), old));
    v(1, new) = v(ceil(width)+1, new);

    temp = old;
    old = new; 
    new = temp;    
end

% Reduces array v to contain only solution
v = v([1:ceil(width), old]);
solution = zeros(1, ceil(width)+1);
error = zeros(1, ceil(width)+1);

% Because the true solution is known, we may plot it. This portion appends
% the true solution to v so that it may be plotted and compared. It would
% not normally be included in this algorithm.
% Additionally, makes an array to use as x-axis values
xVals = zeros(ceil(width)+1, 1);
for i = 1: ceil(width) + 1
    solution(i) = u_0((x_min + (i-1)*h) - a*(t_max));
    xVals(i) = x_min + (i-1) * h;
    error(i) = abs(v(i) - solution(i));
end

% Flipping the vector to make it easy to plot
v = transpose(v);

plot(xVals, v, "--", xVals, solution, "-")
title("Lax-Wendroff")
xlabel("x")
legend({'approximation', 'true'}, Location="northwest")

% Removes last element from error vector to not double count endpoints
error = error(1:ceil(width));

norm(error) / sqrt (length(error))
max(error)

end
