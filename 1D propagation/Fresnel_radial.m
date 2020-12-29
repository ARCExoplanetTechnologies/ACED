% Direct Fresnel integration for the case of radially symmetrical functions
%
% This function requires at least F*4pi samples per aperture, where F = a^2/(lambda*Z)). 
% It works well for longer propagations where F is smaller. For larger F,
% consider using the AS version of this routine


function out = Fresnel_radial(Ein, r, q, lambda, z);
k = 2*pi/lambda;
out = (1/(i*lambda*z))*exp(i*k*(z + q.^2/(2*z))).*hankel_transform_memory(Ein.*exp(i*k*r.^2/(2*z)), r, q/(lambda*z)); 