%% (a)
%%
x1 = -10:0.01:10;
y1 = -10.*cos(3.*x1).^2-(x1-5).^2+250;
figure();
plot(x1,y1);
xlabel('x');
ylabel('f(x)');
title('(a)');
%%

%% (b)
%%
X2 = -3:0.1:3;
Y2 = -3:0.1:3;
[x2,y2] = meshgrid(X2,Y2);
z2 = 3.*(1-x2).*(1-x2).*exp(-x2.^2-(y2+1).^2) - 10.*(x2./5-x2.^3-y2.^5).*exp(-x2.^2-y2.^2) - 1/3.*exp(-(x2+1).^2-y2.^2);
figure();
surf(x2,y2,z2);
xlabel('x');
ylabel('y');
zlabel('f(x,y)');
title('(b)');
%%

%% (c)
%%
X3 = -1:0.02:1;
Y3 = -1:0.02:1;
[x3,y3] = meshgrid(X3,Y3);
z3 = x3.^2 + y3.^2 - 0.5.*cos(pi.*x3) - 0.5*cos(2.*pi.*y3) + 1;
figure();
surf(x3,y3,z3);
xlabel('x');
ylabel('y');
zlabel('f(x,y)');
title('(c)');
%%

%% (d)
%%
X4 = -1:0.02:1;
Y4 = -1:0.02:1;
[x4,y4] = meshgrid(X4,Y4);
z4 = x4.^2 + y4.^2 - 0.7.*cos(2.*pi.*x4).*cos(3.*pi.*y4) + 0.7;
figure();
surf(x4,y4,z4);
xlabel('x');
ylabel('y');
zlabel('f(x,y)');
title('(d)');
%%

%% (e)
%%
X5 = -1:0.02:1;
Y5 = -1:0.02:1;
[x5,y5] = meshgrid(X5,Y5);
z5 = x5.^2 + 2.*y5.^2 - 0.3.*cos(4.*pi.*x5) - 0.3.*cos(5.*pi.*y5) + 0.6;
figure();
surf(x5,y5,z5);
xlabel('x');
ylabel('y');
zlabel('f(x,y)');
title('(e)');
%%