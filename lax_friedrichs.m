function [v] = lax_friedrichs(h, lambda)
% A numerical solution to the one-way wave equation. This numerical
% solution uses forward steps in time, and central steps in space.

syms x;

% Constants. These could be specified by user input but I just hard coded
% them in
x_min = -1;
x_max = 3;
t_min = 0;
t_max = 2.4;

% Technically, this is specified by the user through inputting h and lambda
k = lambda * h;

% Initial data
u_0(x) = piecewise(abs(x) <= 0.5, (cos(pi * x))^2, 0);

% Create array to store function values and fill lowest time level with
% initial data
width = (x_max - x_min) / h;
height = (t_max - t_min) / k;
v = zeros(ceil(width)+1, 2);
for i = 1 : ceil(width)+1
    v(i,1) = u_0(x_min + h * (i-1));
end

% Note that the matrix v has position (x) as the rows, and time as the
% columns, i.e. (x, t)

old = 1;
new = 2;

% Compute next layer of time values
for j = 1 : ceil(height)
    % Compute all space values at given time before moving on
    for i = 2 : ceil(width)
        % Formula for Lax-Friedrichs
        v(i, new) = -k * (v(i+1, old) - v(i-1,old))/(2*h) + (1/2)*(v(i+1, old) + v(i-1, old));
    end
    % Right boundary condition
    v(ceil(width)+1, new) = v(ceil(width), new);

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
    solution(i) = u_0((x_min + (i-1)*h) - (t_max));
    xVals(i) = x_min + (i-1) * h;
    error(i) = abs(v(i) - solution(i));
end

% Flipping the vector to make it easy to plot
v = transpose(v);

plot(xVals, v, "--", xVals, solution, "-", xVals, error, "-.")
title("Lax-Friedrichs")
xlabel("x")
legend({'approximation', 'true', 'error'}, Location="northwest")

norm(error)

end