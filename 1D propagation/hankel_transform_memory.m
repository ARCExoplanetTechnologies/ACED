% This function evaluates the hankel transform of Fin
% Fin is assumed to be evaluated at values r
% r is assumed to be a row vector, uniform array of values starting from 0
% output is assumed to be evaluated at q, also a row vector

function out = hankel_transform_memory(Fin, r, q);

N = length(r);
% dr(1) = (r(2) - r(1))/2;
% dr(2:N-1) = (r(3:N) - r(1:N-2))/2;
% dr(N) = r(N) - r(N-1);
dr = r(2) - r(1);

out = zeros(1, length(q));

rdr = r*dr;
%rdr(1) = dr^2/8;
%int_rule = ones(1,N); % rectangle rule
int_rule = ones(1,N); int_rule(1) = 1/2; int_rule(N) = 1/2; % trapezoid rule
%int_rule = (mod((0:(N-1)),2)+1)*2/3; int_rule(1) = 1/3; int_rule(N) = 1/3; % Simpson's rule, requires N to be even;

for i = 1:length(q)
    out(i) = 2*pi*besselj(0, 2*pi*q(i)*r)*((Fin.*r.*int_rule*dr).'); % rectangle rule integration
end