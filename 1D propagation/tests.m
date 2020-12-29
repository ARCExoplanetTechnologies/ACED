clear all;
N = 512;
M = N;
lambda = 632e-9;
f = 1;
D = 16e-3;

flD = f*lambda/D;

dr = D/N;
r = (0:(N-1))*dr;
q = (0:(M-1))/M*10*flD;

Ein = r < (D/2);

% tic; Fout = hankel_transform(Fin, r, q); toc;
% bellerophon times:
% N,M      time
% 1024     6s
% 2048     19s
% 4096     72s
% 8192     out of memory

% tic; Fout = hankel_transform_memory(Fin, r, q); toc;
% bellerophon times:
% N,M      time
% 1024     3s
% 2048     18s
% 4096     73s
% 8192     304s

% *** Fresnel radial and radial angular spectrum test ***

z = 1;
tic; Eout1 = Fresnel_radial(Ein, r, r, lambda, z);toc;
tic; Eout2 = Fresnel_radial_AS(Ein, r, r, lambda, z, 1/2, 2); toc;

figure(1)

plot(r, abs(Ein).^2,'b');hold on;
plot(r, abs(Eout1).^2, 'g');
plot(r, abs(Eout2).^2, 'r'); hold off;
axis([0 D 0 2]);

Fresnel_number = (D/2)^2/(lambda*z)
%Minimum number of samples required
Min_samples = Fresnel_number*pi
Validity_of_Fresnel = z^3/(pi/(lambda) * D^4)



% 
% % LENS TEST
% % Eout_far = Fraunhofer_radial(Fin, r, q, lambda, z);
% E1 = Fresnel_radial(Ein, r, r, lambda, 2*f);
% E2 = E1.*lens(f, lambda, 2*D, r);
% figure(1)
% plot(r, abs(E1).^2);
% 
% Eout_near = Fresnel_radial(E2, r, q, lambda, 1);
% 
% Fresnel_number = (D/2)^2/(lambda*f)
% %Minimum number of samples required
% Fresnel_number*pi
% 
% figure(2);
% %plot(r, abs(Fin).^2,'b'); hold on;
% %plot(q, abs(Eout_far).^2,'g');
% plot(r, abs(Eout_near).^2,'r');
% axis([0 0.001 0 2]);