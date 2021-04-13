function [y,g] = phi(lambda, t, x)
    z = sqrt(lambda^2 * x^2 + t^2);
    y = lambda * x + z;
    g = lambda + (x*lambda^2) / z; 
end