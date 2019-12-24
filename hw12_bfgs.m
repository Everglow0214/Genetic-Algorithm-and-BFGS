clear all

x0 = [-1,-1]';
% x0 = [1,-0.75]';

rho = 0.55;
sigma = 0.4;
epsilon = 1e-5;

k = 0;
iteration = 1000;
n = length(x0);
Bk = eye(n);
while(k < iteration)
    %gk = feval(gfun, x0);
    gk = gfun(x0);
    if (norm(gk) < epsilon)
        break;
    end
    dk = -Bk\gk;
    m = 0;
    mk = 0;
    while(m < 20)
        %newf = feval(fun, xo+dk*rho^m);
        %oldf = feval(fun, x0);
        newf = fun(x0+dk*rho^m);
        oldf = fun(x0);
        if (newf < oldf+sigma*rho^m*(gk'*dk))
            mk = m;
            break;
        end
        m = m+1;
    end
    x = x0+dk*rho^m;
    sk = x-x0;
    %yk = feval(gfun,x)-gk;
    yk = gfun(x)-gk;
    if (yk'*sk>0)
        Bk = Bk- (Bk*sk*sk'*Bk)/(sk'*Bk*sk) + (yk*yk')/(yk'*yk);
    end
    k = k+1;
    x0 = x;
end

x0
fun(x0)

function f=fun(x)
    f = x(1).^2 + 2.*x(2).^2 - 0.3.*cos(4.*pi.*x(1)) - 0.3.*cos(5.*pi.*x(2)) + 0.6;
end

function gf=gfun(x)
    gf = [2*x(1)+1.2*pi*sin(4.*pi.*x(1)),4*x(2)+1.5*pi*sin(5.*pi.*x(2))]';
end