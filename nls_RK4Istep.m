function uRight = nls_RK4Istep(uLeft,gLL,gLeft,gRight,gRR,lfac,nt,dz);

% Integrates u by one step according to forced NLS eqn, i.e.
%
% du/dz = i*lfac*d^2u/dt^2 + i*|u|^2*u + g
%
% over single step dz with t\in (-pi,pi)
%
% Using integrating-factor RK4I method (pseudospectral)
%
% Uses 4th-order interpolation to get g at half-steps

k = [0:nt/2-1 -nt/2:-1]';
E = exp(-0.5i*dz*lfac*k.^2);
E2 = exp(-1i*dz*lfac*k.^2);

uftLeft = fft(uLeft);
gftLeft = fft(gLeft);
gftRight = fft(gRight);
gftLL = fft(gLL);
gftRR = fft(gRR);

uCur = uLeft;
gftCur = gftLeft;
a = dz*1i*abs(uCur).^2.*uCur;
a = fft(a) + dz*gftCur;

uCur = ifft(E.*(uftLeft + a/2));
gftCur = 9/16*(gftLeft + gftRight) - 1/16*(gftLL+gftRR);
b = dz*1i*abs(uCur).^2.*uCur;
b = fft(b) + dz*gftCur;

uCur = ifft(E.*uftLeft + b/2);
c = dz*1i*abs(uCur).^2.*uCur;
c = fft(c) + dz*gftCur;

uCur = ifft(E2.*uftLeft + E.*c);
gftCur = gftRight;
d = dz*1i*abs(uCur).^2.*uCur;
d = fft(d) + dz*gftCur;

uftRight = E2.*uftLeft + (E2.*a + 2*E.*(b+c) + d)/6;
uRight = ifft(uftRight);

