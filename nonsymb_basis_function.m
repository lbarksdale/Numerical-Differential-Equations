function p = nonsymb_basis_function(h, i, x_min, xvals)
% Creates an array of basis function values given input of h (mesh size),
% i (index of basis function), x_min (the left side of the domain), and
% xvals (an array of x values)

p = zeros(length(xvals), 1);

for k = 1:length(xvals)
    if(x_min + (i-1) * h <= xvals(k) && xvals(k) <= x_min + i*h)
        p(k) = 1/h * (xvals(k) - (i-1)*h);
    elseif(x_min + (i)*h <= xvals(k) && xvals(k) <= x_min + (i+1)*h)
        p(k) = -1/h * (xvals(k) - (i+1)*h);
    else
        p(k) = 0;
    end
end


end
