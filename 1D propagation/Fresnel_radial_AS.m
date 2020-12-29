% This function evaluates the radially symmetrical Fresnel propagation
% using the Angular spectrum method
%
% This routine suffers from artifacts due to numerical integration (at least if the sampling theorem based on bessel-function-zeros is not used), 
% but gives a "better" answer than direct Fresnel integration for the case of insufficient samples. 
% ASF effectively gives a low-passed-filtered version of the true result for the case that there are insufficient samples 

% 2 parameters to control the numerical effects -- a LPF factor (which is 1 for "Nyquist"  frequency), and an oversampling factor of the Hankel space. 
% The product of these is the ratio of the number of samples in the hankel space to the number of samples in the physical space.
% In other words, higher LPF and Oversampling values slow down computations


function out = Fresnel_radial_AS(Ein, r, q, lambda, z, LPF, Oversampling);
k = 2*pi/lambda;
N = length(r);

dr = r(2) - r(1); % Assumes uniform sampling of r
M = N*Oversampling*LPF; % M = N would be the FFT sampling rate, but hankel obeys a different sampling theorem (M needs to be >> N for more accuracy)
dnu = LPF*(0.5/dr)/M; % go up to the "Nyquist" frequency
nu = (0:(M-1))*dnu;

H = exp(i*k*z * (sqrt( 1 - (lambda*nu).^2)));

out = hankel_transform_memory( hankel_transform_memory(Ein, r, nu).*H, nu, q);