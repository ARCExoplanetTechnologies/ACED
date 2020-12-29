% Simulations of free-space fresnel propagation from a square aperture as seen on
% page 86 of Goodman 3rd ed.
%
% Here are a few demonstration cases to try, going from large Fresnel
% numbers (close in) to small Fresnel numbers (far out). All these cases
% assume a square aperture 0.5mm on the side (a = 0.00025).
%
% Case 1: Nf = 100, N = 512, Din = 1e-3 (or anything higher)
% This corresponds to about a 1mm travel distance
% AS is the best method to use here, and as always it assumes an infinite
% periodic plane of apertures, so some zeropadding is necessary (Try
% reducing Din until there is no zeropadding to see what happens)
%
% The reason the direct method doesn't work for this case is because the quadratic phase
% factor (see FresnelPropagate.m code) is not sampled well enough to have at least 1 sample per radian of
% phase. To prevent this, we need at least 4pi*Nf = 1256 samples in the input sampling
% grid. For example, try the following case:
% 
% Case 2: Nf = 100, N = 2048, Din = 1e-3
% Now the direct method works. However, note how much faster the AS is.
% Also, with the AS we can make N almost as small as we want to improve speed further and it will still give the
% right answer. One advantage of the direct method though is that it doesn't assume the
% infinitely periodic grid of apertures like AS does and so doesn't require
% any zeropadding:
%
% Case 3: Nf = 100, N = 1024, Din = 0.5e-3 (no zeropadding). Note the
% failure of the AS method now, while the direct method works well without
% zeropadding.
%
% In all these cases, the Fraunhofer method fails misearably. Or rather,
% what it effectively calculates is the propagation of our square aperture, but starting with a curved wavefront, 
% i.e. it implicitly assumes that there is a lens in the plane of our aperture with a focal length
% z (i.e. ~ 1mm). Such a lens will focus our 0.5mm square beam 1mm away into a very small 2D sinc pattern with
% a very high peak brightness, which is exactly what we get if you look at
% figure 5. For propagating to foci of lenses, Fraunhofer succeeds admirably, but
% not for propagations to intermediate planes.
%
% Case 4: Nf = 10, N = 512, Din = 1e-3. This corresponds to about a 1cm propagation distance. Both the AS and the direct
% method give similar results. Nf is still comfortably large for AS and 512
% samples is way more than enough to satisfy N > 4pi*Nf for the direct
% method. Fraunhofer still fails of course because we're not in far field.
% Note how much faster AS is than the other two.
%
% Case 5: Nf = 1, N = 512, Din = 1e-3. This corresponds to about a 10 cm propagation distance. Now the AS method starts visibly
% failing, and the direct method is great. In order to bring the accuracy
% of the AS method back, the input needs to be zeropadded more (e.g. Din =
% 8e-3, N = 4096). However, with this much zeropadding, its speed advantages vanish (as compared to the non-zeropadded direct method).
%
% Case 6: NF = 0.1, N = 512, Din = 1e-3. This corresponds to about a 1m
% propagation distance, well after the beam starts diverging, so the Fraunhofer
% approximation becomes accurate, matching the direct method and showing the
% top of the 2D sinc pattern (equivalent to Airy pattern, but for square apertures). 
% For such large distances, there is
% practically no difference between a 1m free-space propagation or a
% propagation after passing through a lens with a 1m focal length (which
% looks almost like a flat window to a 0.5mm beam), so the Fraunhofer
% approximation applies to both cases.
% 
% However, the AS method now fails unless the zeropadding is really large, for example 
%
% Case 7: Nf = 0.1, N = 512, Din = 2e-2. In this case all three methods
% agree and generate the 2D sinc pattern characteristic of the rectangular
% aperture. For this case, zeropadding is actually nice because the output
% field is so much larger than the input field, but for the direct methods
% and Fraunhofer, we can have separate sampling grids at the input and the
% output, thus preserving the large field of view at the output while not zeropadding at the input
% (for AS this is impossible because they are tied together):
%
% Case 8: Nf = 0.1, N = 512, Din = 0.5e-3, but change Dout on line 80 to 2e-2. 

clear all; clc;

N = 512;
% Fresnel number
Nf = 10;
% input plane definition
Din = 1e-3;

x = ((1:N)-N/2)/N*Din;
y = x;

[xx yy] = meshgrid(x,y);
rr = sqrt(xx.^2 + yy.^2);
ttheta = atan2(yy,xx);

% Defining the input field
a = 0.00025; %radius of the square aperture
Ein = abs(xx) < a & abs(yy) < a; % square aperture
Einx = Ein(N/2,:);
%Ein = rr < a; % circular aperture

lambda = 650e-9;
% Propagation distance
z = a^2/lambda/Nf

% output plane definition (ignored by the AS method, where Dout = Din always)
Dout = Din;
xi = ((1:N)-N/2)/N*Dout;
eta = xi;

[xxi eeta] = meshgrid(xi,eta);
rrho = sqrt(xxi.^2 + eeta.^2);
pphi = atan2(xxi,eeta);

% Fresnel, "Direct" method
fprintf('Computing with Direct method...\n');
tic
EoutD = FresnelPropagate(x,y,Ein,xi,eta,z,lambda);
EoutDx = EoutD(N/2,:);
toc

% Fresnel, Angular spectrum method
fprintf('Computing with AS method...\n');
tic
EoutAS = FresnelPropagateAS(Ein,lambda,Din/2,z);
EoutASx = EoutAS(N/2,:);
toc

% Fraunhofer
fprintf('Computing Fraunhofer...\n');
tic
EoutF = zoomFFT_realunits(x,y,Ein,xi,eta,z,lambda);
EoutFx = EoutF(N/2,:);
toc

figure(1)
imagesc(x/1e-3,y/1e-3,abs(Ein).^2); axis image;
axis([-3 3 -3 3]);
title 'Input Field';

figure(2)
imagesc(x/1e-3,y/1e-3,abs(EoutAS).^2); axis image; colorbar;
axis([-3 3 -3 3]);
str = sprintf('Output Field Nf= %f',Nf);
title('Angular Spectrum Method');

figure(3)
imagesc(xi/1e-3,eta/1e-3,abs(EoutD).^2); axis image; colorbar;
axis([-3 3 -3 3]);
str = sprintf('Output Field Nf= %f',Nf);
title('Direct method');

figure(4)
imagesc(xi/1e-3,eta/1e-3,abs(EoutF).^2); axis image; colorbar;
axis([-3 3 -3 3]);
str = sprintf('Output Field Nf= %f',Nf);
title('Fraunhofer');

figure(5)
semilogy(x/1e-3, abs(Einx).^2+1e-10,'k'); hold on;
semilogy(xi/1e-3, abs(EoutDx).^2); hold on;
semilogy(x/1e-3, abs(EoutASx).^2,'g'); 
semilogy(xi/1e-3, abs(EoutFx).^2,'r'); hold off;
str = sprintf('Output Field Nf= %f',Nf);
title(str);
legend('Input Field', 'Fresnel, direct','Fresnel, AS', 'Fraunhofer');