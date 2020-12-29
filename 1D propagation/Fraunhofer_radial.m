% This function evaluates the hankel transform of Fin
% Fin is assumed to be evaluated at values r
% r is assumed to be a row vector, uniform array of values starting from 0
% output is assumed to be evaluated at q, also a row vector

function out = Fraunhofer_radial(Ein, r, q, lambda, z);
k = 2*pi/lambda;
out = (1/(i*lambda*z))*exp(i*k*(z + q.^2/(2*z))).*hankel_transform_memory(Ein, r, q/(lambda*z)); 