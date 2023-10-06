%Script to test Galerkin method
syms x;
f(x) = sym(0);
n = 20;
x_min = 0;
x_max = 1;
h = (x_max - x_min) / n;
xvals = linspace(x_min, x_max, n+1);
solxvals = x_min : h/10: x_max;
u = galerkin_finite_element(n, x_min, x_max, f(x));

dense_approx = interp1(xvals, u, solxvals);

% Create exact solution plot

c2 = 1/(sqrt(2) * (besselk(0,2) - besselk(0,2*sqrt(2)) * besseli(0,2) / besseli(0, 2*sqrt(2))));
c1 = -c2 * besselk(0, 2*sqrt(2)) / besseli(0,2*sqrt(2));

solyvals = zeros(length(solxvals),1);

for i = 1:length(solyvals)
    solyvals(i) = sqrt(2) * (c1 * besseli(0, 2 * sqrt(solxvals(i) + 1)) + c2 * besselk(0, 2 * sqrt(solxvals(i) + 1))); 
end

plot(solxvals, dense_approx, solxvals, solyvals, "--")

error = abs(transpose(dense_approx) - solyvals);
format longg
double(max(error))
legend('approximate', 'true')
figure
plot(solxvals, error)

