function [y,g] = theta(t, x)
    y = sqrt(sum(x.*x) + t^2);
    g = (1/y) * x;
end