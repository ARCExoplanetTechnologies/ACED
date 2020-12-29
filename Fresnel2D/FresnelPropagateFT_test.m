N = 1024;
lambda = 650e-9;
z = 2;

% input plane definition
a = 2e-2/2;

x = ((1:N)-N/2)/N*2*a;
y = x;

[xx yy] = meshgrid(x,y);
rr = sqrt(xx.^2 + yy.^2);
ttheta = atan2(yy,xx);

Ein = rr < 5e-4;

Eout = FresnelPropagateFT(Ein, lambda, a, z);

figure(1)
imagesc(x,y,abs(Ein).^2); axis image;

figure(2)
imagesc(x,y,log((abs(Eout)).^2),[-8 2]); axis image;colorbar;