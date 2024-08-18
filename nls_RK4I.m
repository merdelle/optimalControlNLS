function [uRight,uftRight] = nls_RK4Istep(uLeft,gLeft,gRight,lfac,nt,dz);

% Integrates u by one step according to forced NLS eqn, i.e.
%
% du/dz = i*lfac*d^2/dt^2 + i*|u|^2*u + g
%
% over single step dz with t\in (-pi,pi)

dt = .1/(2*N);
x = (2*pi/N)*(-N/2:N/2-1)';
% Set initial condition
lam = 10;
u = sqrt(2) * lam * sech(lam * x);
v = fft(u);

% Construct integration operators
k = [0:N/2-1 0 -N/2+1:-1]';
E = exp(-.5i * dt * k.^2);
E2 = E.^2;
g = 1i * dt;

% Set simulation parameters
tmax = 5; nplt = floor((tmax/50)/dt); nmax = round(tmax/dt);
udata = u; tdata = 0;

% Integrate using the Integrating Factor method
for n = 1:nmax
  t = n*dt;
  a = g.*fft(abs( ifft(     v   ) ).^2 .* ifft(v));
  b = g.*fft(abs( ifft(E.*(v+a/2)) ).^2 .* ifft(E.*(v+a/2)) );
  c = g.*fft(abs( ifft(E.*v + b/2) ).^2 .* ifft(E.*v + b/2) );
  d = g.*fft(abs( ifft(E2.*v+E.*c) ).^2 .* ifft(E2.*v+E.*c));
  v = E2.*v + (E2.*a + 2*E.*(b+c) + d)/6;
  if mod(n,nplt) == 0
    u = ifft(v);
    udata = [udata u]; tdata = [tdata t];
  end
end
