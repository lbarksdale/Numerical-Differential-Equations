function [area]  = simpson(vals, step)
% Simple function to numerically integrate using composite Simpson's method. vals is
% an array of y-values, and step is the spacing in x
% Note that every other x-value is considered a "midpoint" when using
% Simpson's

area = 0;

for i = 1:2:length(vals) - 2
    area = area + step / 3 * (vals(i) + 4 * vals(i+1) + vals(i+2));
end

end