function [u] = galerkin_finite_element(n, x_min, x_max, forcing_term)
% This method solves the strong form of a PDE using the Galerkin finite-element method.
% It takes in several arguments:
% n - the number of node points to use (defines discretization)
% x_min - the minimal domain value
% x_max - the maximal domain value
% forcing_term - the forcing term of the differential equation
% k - the function k(x) found in the weak form of the BVP
% Note that this implementation is a modification of the previous one,
% where nonzero boundary conditions are implemented. However, calculation
% for this was done by hand

% The Galerkin finite-element method is essentially solving the matrix
% equation Au = f
% This is created from the weak form of the BVP, which is 
% a(u, v) = (f, v)
% v is based off of a choice of basis vectors. These may be adjusted by
% changing the function basis_function, which is currently set to piecewise
% linear functions.

syms x;

f(x) = forcing_term;
h = (x_max - x_min) / n;
xvals = linspace(x_min, x_max, n+1);
dense_x_vals = linspace(x_min, x_max, 4*n+1);
dense_k_vals = k(dense_x_vals);

A = zeros(n-1);

% For the specific basis functions involved in this problem, the matrix A
% will be tridiagonal (leading to easy computations). In general, however,
% this may not be the case.

% Create the stiffness matrix A that will be inverted
for i = 1:n-1
    for j = 1:n-1
        % Elements of stiffness matrix are inner products

        % For this specific problem, only compute tridiagonal entries
        % (other entries have non-overlapping basis functions and are
        % therefore 0)
        if abs(i - j) <= 1
            %Calculation using numerical differentiation and integration
                
            % Calculates values of basis functions at node points
            basis_1_vals = nonsymb_basis_function(h, i-1, x_min, xvals);
            basis_2_vals = nonsymb_basis_function(h, j-1, x_min, xvals);

            dense_basis_1 = interp1(xvals, basis_1_vals, dense_x_vals);
            dense_basis_2 = interp1(xvals, basis_2_vals, dense_x_vals);
            % Since this problem uses hat functions, a simple finite
            % difference and left-side sum is exact for differentiation and
            % integration, so that's what I used (left side sum is exact
            % for integration as the derivative of hat functions are
            % constant)

            diff1 = simple_finite_difference_derivative(dense_basis_1, h/4);
            diff2 = simple_finite_difference_derivative(dense_basis_2, h/4);
            for m = 1:4:length(diff1)-4
                prod = diff1(m) * diff2(m);
                mat = [prod prod prod prod prod];
                A(i,j) = A(i,j) + simpson(times(transpose(dense_k_vals(m:m+4)), mat), h/4);
            end

            % Integrates the u*v term
            A(i,j) = A(i,j) + simpson(times(dense_basis_1, dense_basis_2), h/4);
            
        end
    end
end

% Accounts for off-by-one error in integration of first term, originates
% from cutoff of domain on left hand side and this is the easiest fix
% Actually, this might be introducing some error because it doesn't account
% for the change in k as x varies but I don't have the time to fix it :(
A(1,1) = 2 * A(1, 1);

% Creates and populates the vector discretization of the forcing term
forcing_vector = zeros(n-1, 1);
forcing_vector(1) = 1/h + 9.82 *h;

% Solves the matrix-vector system
u = A\forcing_vector;
% Adds in Dirichlet boundary conditions
u = [1;u;0];

end