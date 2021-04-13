clear all
ntest = 100;
t = rand(ntest, 1);
lambda = rand(ntest, 1);
x = -1 + 2 * rand(ntest, 1);

for i = 1:ntest
    [~,gx] = phi(lambda(i), t(i), x(i)); 
    gn = num_diff(@(x)phi(lambda(i), t(i), x), x(i));
    
    if abs(gn - gx) > 1E-3
        error('num_diff is different')
    end
end