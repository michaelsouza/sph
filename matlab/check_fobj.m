function check_fobj()
test_A();
test_B();
end

function test_A()
[G,x] = create_grid_instance(5,5,0.1,0);
x = x(:);
[fx,gx] = fobj(G,x);
if abs(fx) > 1E-8
    error('The fx should be zero.')
end

if norm(gx) > 1E-8
    error('The gx should be null.')
end

y = x;
y(end) = 1;
[~,gy] = fobj(G,y);
gn = num_diff(@(y)fobj(G,y),y);
if norm(gy - gn) > 1E-3
    error('The num_diff is different')
end
end

function test_B()
[G,x] = create_grid_instance(20,5,0.1,0);
x = x(:);
ntest = 100;
for i = 1:ntest
    y = rand(size(x));
    [~,gy] = fobj(G,y);
    gn = num_diff(@(y)fobj(G,y),y);
    if norm(gy - gn) > 1E-3
        error('The num_diff is different')
    end
end
end