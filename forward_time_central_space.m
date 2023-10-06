function [v] = forward_time_central_space(h, max_time_layers)
% A numerical solution to the one-way wave equation. This numerical
% solution uses forward steps in time, and central steps in space.

syms x;

% Allows user to specify time layers value of -1 if they want to simply
% compute to the end of the time domain
if max_time_layers == -1
    max_time_layers = intmax;
end

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
% initial data
width = (x_max - x_min) / h;
height = (t_max - t_min) / k;
v = zeros(ceil(width)+1, 2);
for i = 1 : ceil(width)
    v(i,1) = u_0(x_min + h * (i-1));
end

% Note that the matrix v has position (x) as the rows, and time as the
% columns, i.e. (x, t)

old = 1;
new = 2;

useless = false;

% Compute next layer of time values
for j = 1 : ceil(height)-1
    % Allows for user to specify the number of time layers that should be
    % computed
    if j > max_time_layers
        break
    end
    % Compute all space values at given time before moving on
    for i = 2 : ceil(width)-1
        % Formula for forward time, central space
        v(i, new) = -k * (v(i+1, old) - v(i-1,old))/(2*h) + v(i, old);
        if v(i, new) > 5
            useless = true;
        end
    end
    % Right boundary condition
    v(ceil(width), new) = v(ceil(width) - 1, new);

    temp = old;
    old = new;
    new = temp;
end

% Reduces array v to contain only solution
v = v([1:ceil(width), old]);
solution = zeros(1,ceil(width)+1);
error = zeros(1,ceil(width)+1);

% Because the true solution is known, we may plot it. This portion appends
% the true solution to v so that it may be plotted and compared. It would
% not normally be included in this algorithm.
% Additionally, makes an array to use as x-axis values
xVals = zeros(ceil(width)+1, 1);
for i = 1: ceil(width) + 1
    solution(i) = u_0((x_min + (i-1)*h) - (t_min + max_time_layers * k));
    error(i) = abs(v(i) - solution(i));
    xVals(i) = x_min + (i-1) * h;
end

% Flipping the vector to make it easy to plot
v = transpose(v);

plot(xVals, v, "--", xVals, solution, "-", xVals, error, "-.")
if(useless)
    title("Forward-Time Central-Space (USELESS)")
else
    title("Forward-Time Central-Space")
end
xlabel("x")
legend({'approximation', 'true', 'error'}, Location="northwest")

end