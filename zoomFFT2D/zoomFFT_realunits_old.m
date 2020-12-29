function out = zoomFFT_realunits(f,lambda,a,Npup,u,Nimg,Ein)

% Fraunhover integral (f-f lens case, with e^(i*k*z)/i factor dropped)
% 
% Computation is accomplished via a zoomFFT (chirpZ).
%
% Ein and the output are both assumed to be defined on a square symmetric
% about the origin
%
% a = diameter of aperture in meters
% 2*Npup = number of points in pupil plane
% u = size of image plane in meters(CCD)
% 2*Nimg+1 = number of points in image plane (pixels in CCD)

dxi = u/2/Nimg;
deta = u/2/Nimg;
dx = a/2/Npup;
dy = a/2/Npup;

xs = [-Npup:Npup-1]'*dx + dx/2;
ys = [-Npup:Npup-1]'*dy + dy/2;

xis = [-Nimg:Nimg]'*dxi;
etas = [-Nimg:Nimg]'*deta;

A = exp(2*pi*i * dx * xis(1)/(lambda*f));
W = exp(-2*pi*i * dx * dxi/(lambda*f));

% sorry about this line, I was trying to optimize...
out = czt(czt(Ein, 2*Nimg+1, W, A).',2*Nimg+1, W, A).'.* exp(-2*pi*i * (xs(1)*ones(2*Nimg+1,1)*xis' + ys(1)*etas*ones(1,2*Nimg+1))/(lambda*f)) * dx*dy /(lambda*f);

% Below are various tests

% A = exp(-pi*i * a * u /(2*Npup*lambda*f));
% W = exp(-pi*i * a * u /(2*Npup*Nimg*lambda*f));
 
% out = czt(czt(Ein, 2*Nimg, W, A).', 2*Nimg, W, A).'.* exp(pi* i * (Npup - 1/2) * a * u * (ones(2*Nimg,1)*[-Nimg:Nimg-1] + [-Nimg:Nimg-1]'*ones(1,2*Nimg)) /(2*Nimg*Npup*lambda*f)) * a^2 /(4*Npup^2*lambda*f);
% out = czt(czt(Ein, 2*Nimg, W, A).', 2*Nimg, W, A).'.* exp(i*(ones(2*Nimg,1)*[-Nimg:Nimg-1]));