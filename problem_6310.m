function [v] = problem_6310(h, lambda, mu)
% A numerical solution to the one-way wave equation. This numerical
% solution uses forward steps in time, and backwards steps in space.

% This program is adapted from an earlier assignment, but modified to
% better demonstrate order of accuracy.

syms x;

if lambda == 0
    k = mu * h^2;
    lambda = k / h;
else
    k = lambda * h;
    mu = k / h^2;
end

% Constants
x_min = -1;
x_max = 1;
t_min = 0;
t_max = 0.5;
b = 1;

% Initial data
u_0(x) = piecewise(abs(x)<1/2, 1, abs(x)==1/2, 1/2, 0);

% Create array to store function values and fill lowest time level with
% initial data. Array consists of two levels, an old and a new
width = (x_max - x_min) / h;
height = (t_max - t_min) / k;
v = zeros(ceil(width)+1, 2);
for i = 1 : ceil(width)+1
    v(i,1) = u_0(x_min + h * (i-1));
end

% Creates matrix of coefficients of linear system that we solve (it's an
% implicit scheme, so we have to solve a linear system to get the solution
% at the next time value)
coefficients = zeros(width - 1);
r_coeffs = zeros(width - 1);
r_coeffs(1,1) = 1 - b * mu;
r_coeffs(1,2) = b * mu / 2;
coefficients(1,1) = 1 + b * mu;
coefficients(1,2) = -b * mu / 2;

for i = 2:width-2
    coefficients(i, i-1) = -b * mu / 2;
    coefficients(i, i) = 1 + b * mu;
    coefficients(i, i+1) = -b * mu / 2;
    r_coeffs(i, i-1) = b * mu / 2;
    r_coeffs(i,i) = 1 - b * mu;
    r_coeffs(i, i+1) = b * mu / 2;
end

coefficients(width - 1, width - 2) = -b * mu / 2;
coefficients(width - 1, width - 1) = 1 + b * mu;
r_coeffs(width - 1, width - 2) = b * mu / 2;
r_coeffs(width - 1, width - 1) = 1 - b * mu;

% Old and new track which layer of the array v holds old data vs new data
old = 1;
new = 2;

% Keep track of computation time
tic;
% Compute next layer of time values until desired time level is reached
for j = 1 : ceil(height)
    % Compute all space values at given time before moving on
    % Writes new data over old (unneeded) data

    % Computes forcing term (in this problem forcing term is 0, but it can
    % easily be added here if it is not)
    f = zeros(width -1, 1);
    f_newer = zeros(width - 1, 1);

    % Boundary Data
    v(1, new) = u(t_min + j*k, x_min, 4);
    v(ceil(width)+1, new) = u(t_min + j*k, x_max, 4);

    % Additional coefficients on start and end term
    edges = zeros(ceil(width) - 1, 1);
    edges(1) = b * mu / 2 * (v(1, old) + v(1, new));
    edges(width - 1) = b * mu / 2 * (v(width + 1, old) + v(width + 1, new));

    % Creates the linear system to solve
    RHS = r_coeffs*v(2:end-1, old);
    RHS = RHS + k * (f + f_newer) / 2 + edges;
    new_data = linsolve(coefficients, RHS);

    v(2:end - 1, new) = new_data;

    temp = old;
    old = new;
    new = temp;    
end

% Ends recorded time of computation
time = toc;

% Reduces array v to contain only solution
v = v([1:ceil(width), old]);
solution = zeros(1, ceil(width)+1);
error = zeros(1, ceil(width)+1);

% Because the true solution is known, we may plot it. This portion appends
% the true solution to v so that it may be plotted and compared. It would
% not normally be included in this algorithm.
xVals = zeros(ceil(width)+1, 1);
for i = 1: ceil(width) + 1
    solution(i) = u(t_max, x_min + (i-1)*h, 3);
    xVals(i) = x_min + (i-1) * h;   % X-axis values
    error(i) = abs(v(i) - solution(i));
end

% Flipping the vector to make it easy to plot
v = transpose(v);

plot(xVals, v, "--", xVals, solution, "-")
title("Crank-Nicolson Solution to Heat Equation")
xlabel("x")
ylabel("u")
legend({'approximation', 'true'}, Location="northwest")

disp(['Computation time: ', num2str(time)])
disp(['L2 Error: ', num2str(norm(error) / sqrt(length(error)))])
disp(['Supremum Error: ', num2str(max(error))])

end