function gLeft = lnls_RK4Istep(gRight,uLL,uLeft,uRight,uRR,lfac,nt,dz);

% Integrates g by one backward step according to linearized NLS eqn, i.e.
%
% dg/dz = i*lfac*d^2g/dt^2 + 2i*|u|^2*g -i*u^2*conj(g)
%
% over single backward step -dz with t\in (-pi,pi)
%
% Using integrating-factor RK4I method (pseudospectral)
%
% Uses 4th-order interpolation to get u at half-steps

k = [0:nt/2-1 -nt/2:-1]';
E = exp(0.5i*dz*lfac*k.^2);
E2 = exp(1i*dz*lfac*k.^2);

gftRight = fft(gRight);

gCur = gRight;
uCur = uRight;
a = -dz*(2i*abs(uCur).^2.*gCur - 1i*uCur.^2.*conj(gCur));
a = fft(a);

gCur = ifft(E.*(gftRight + a/2));
uCur = 9/16*(uLeft + uRight) - 1/16*(uLL + uRR);
b = -dz*(2i*abs(uCur).^2.*gCur - 1i*uCur.^2.*conj(gCur));
b = fft(b);

gCur = ifft(E.*gftRight + b/2);
c = -dz*(2i*abs(uCur).^2.*gCur - 1i*uCur.^2.*conj(gCur));
c = fft(c);

gCur = ifft(E2.*gftRight + E.*c);
uCur = uLeft;
d = -dz*(2i*abs(uCur).^2.*gCur - 1i*uCur.^2.*conj(gCur));
d = fft(d);

gftLeft = E2.*gftRight + (E2.*a + 2*E.*(b+c) + d)/6;
gLeft = ifft(gftLeft);
