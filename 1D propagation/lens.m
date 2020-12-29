% This function generates a paraxial lens with size D and focal length f
function out = lens(f, lambda, D, r);

out = exp(-i*pi/(lambda*f)*r.^2).*(r<(D/2));