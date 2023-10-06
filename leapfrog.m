function [v] = leapfrog(h)
% A numerical solution to the one-way wave equation. This numerical
% solution uses forward steps in time, and central steps in space.

syms x;

% Constants. These could be specified by user input but I just hard coded
% them in
x_min = -1;
x_max = 3;
t_min = 0;
t_max = 2.4;
lambda = 0.8;
k = lambda * h;

% Initial data
u_0(x) = piecewise(abs(x) <= 0.5, (cos(pi * x))^2, 0);

% Create array to store function values and fill lowest time level with
% initial data
width = (x_max - x_min) / h;
height = (t_max - t_min) / k;
v = zeros(ceil(width)+1, 1);
for i = 1 : ceil(width)+1
    v(i,1) = u_0(x_min + h * (i-1));
end

% Additionally, use forward-time central-space scheme to gather enough
% initial data to use leapfrog scheme
v = [v forward_time_central_space(h, 1) zeros(ceil(width)+1,1)];
% At this point, v consists of all spacial values at two starting times,
% and a row of zeros.

% Note that the matrix v has position (x) as the rows, and time as the
% columns, i.e. (x, t)

old1 = 1;
old2 = 2;
new = 3;

% Compute next time layer until end time is reached
for j = 1 : ceil(height)
    % Compute all space values at given time before moving on
    for i = 2 : ceil(width)-1
        % Formula for Leapfrog
        v(i, new) = -2 * k * (v(i+1, old2) - v(i-1,old2))/(2*h) + v(i, old1);
    end
    % Right boundary condition
    v(ceil(width), new) = v(ceil(width) - 1, new);

    % Rotates through indices in v so that next time layer is written over
    % old unneeded data
    temp = old1;
    old1 = old2;
    old2 = new;
    new = temp;
end

% Reduces array v to contain only solution
v = v([1:ceil(width), old1]);
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
title("Leapfrog")
xlabel("x")
legend({'approximation', 'true', 'error'}, Location="northwest")

norm(error)


end