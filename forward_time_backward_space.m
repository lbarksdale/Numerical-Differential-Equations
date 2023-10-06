function [v] = forward_time_backward_space(h)
% A numerical solution to the one-way wave equation. This numerical
% solution uses forward steps in time, and backwards steps in space.

syms x;

% Constants
x_min = -1;
x_max = 3;
t_min = 0;
t_max = 2.4;
lambda = 0.8;
k = lambda * h;

% Initial data
u_0(x) = piecewise(abs(x) <= 0.5, (cos(pi * x))^2, 0);

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
for j = 1 : ceil(height)
    % Compute all space values at given time before moving on
    % Writes new data over old (unneeded) data
    for i = 2 : ceil(width)+1
        v(i, new) = -k * (v(i, old) - v(i-1,old))/h + v(i, old);
    end

    temp = old;
    old = new; 
    new = temp;    
end

% Reduces array v to contain only solution
v = v([1:ceil(width), old]);
solution = zeros(1, ceil(width)+1);
error = zeros(1, ceil(width)+1);
max_error = 0;

% Because the true solution is known, we may plot it. This portion appends
% the true solution to v so that it may be plotted and compared. It would
% not normally be included in this algorithm.
% Additionally, makes an array to use as x-axis values
xVals = zeros(ceil(width)+1, 1);
for i = 1: ceil(width) + 1
    solution(i) = u_0((x_min + (i-1)*h) - (t_max));
    xVals(i) = x_min + (i-1) * h;
    error(i) = abs(v(i) - solution(i));
    if(error(i) > max_error)
        max_error = error(i);
    end
end

% Flipping the vector to make it easy to plot
v = transpose(v);

plot(xVals, v, "--", xVals, solution, "-", xVals, error, "-.")
title("Forward-Time Backward-Space")
xlabel("x")
legend({'approximation', 'true', 'error'}, Location="northwest")

max_error
norm(error)

end

