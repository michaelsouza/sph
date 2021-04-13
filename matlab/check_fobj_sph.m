function check_fobj_sph()
close all
rng(1)
test_theta(false);
test_phi(true);
test_fobj_sph_A();
test_fobj_sph_B();
end

function test_theta(verbose)
ntest = 100;
for i = 1:ntest
    x = -1 + 2 * rand(3,1);
    t = rand();
    [~,gx] = theta(t,x);
    gn = num_diff(@(y)theta(t,y),x);
    if norm(gx - gn) > 1E-3
        error('The num_diff is different')
    end
end

if verbose
    x = linspace(-1,1);
    y1 = zeros(size(x));
    y2 = y1;
    for i = 1:length(x)
        y1(i) = theta(0.1,x(i));
        y2(i) = theta(0.2,x(i));
    end
    figure; hold on; box on
    plot(x,y1,'LineStyle','-','DisplayName','\theta_{t=0.1}');
    plot(x,y2,'LineStyle','-','DisplayName','\theta_{t=0.2}');
    plot(x,abs(x),'LineStyle','--','DisplayName','|x|');
    legend('show')
end
end

function test_phi(verbose)
ntest = 100;
for i = 1:ntest
    x = rand();
    t = rand();
    lambda = rand();
    [~,gx] = phi(lambda,t,x);
    gn = num_diff(@(y)phi(lambda,t,y),x);
    if norm(gx - gn) > 1E-3
        error('The num_diff is different')
    end
end

if verbose
    x = linspace(-1,1);
    y1 = zeros(size(x));
    y2 = y1;
    y3 = y1;
    for i = 1:length(x)
        y1(i) = phi(0.5,0.1,x(i));
        y2(i) = phi(0.5,0.2,x(i));
        y3(i) = max(x(i),0);
    end
    figure; hold on; box on
    plot(x,y1,'LineStyle','-','DisplayName','\phi_{\lambda=1/2,t=0.1}');
    plot(x,y2,'LineStyle','-','DisplayName','\phi_{\lambda=1/2,t=0.2}');
    plot(x,y3,'LineStyle','--','DisplayName','max(0,x)');
    legend('show','Location','northwest')
end
end


function test_fobj_sph_A()
[G,x] = create_grid_instance(5,5,0.1,0);
x = x(:);
lambda = 1E-8;
t = 1E-8;

[fx,gx] = fobj_sph(lambda,t,G,x);

if norm(gx) > 1E-3 || abs(fx) > 1E-3
    error('The solution is not a critical point')
end
end

function test_fobj_sph_B()
[G,x] = create_grid_instance(20,5,0.1,0);
x = x(:);
ntest = 100;
for i = 1:ntest
    y = -1 + 2 * rand(size(x));
    lambda = rand();
    t = rand();
    [~,gy] = fobj_sph(lambda,t,G,y);
    gn = num_diff(@(y)fobj_sph(lambda,t,G,y),y);
    if norm(gy - gn) > 1E-3
        error('The num_diff is different')
    end
end
end