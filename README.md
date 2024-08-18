# optimalControlNLS
 Matlab script to find optimal control for inhomogeneous nolinear Schrodinger equation
%
% Script to implement optimal control equations on inhomogeneous NLSE
%
% DOB: 21July2018
% Ref: De La Vega et al.
%

% DEBUGGING notes
    % 23Jul18: forward evolution correct (eyeball norm) with
        % soliton evolution, no conv check done yet
% 26Jul18: adopting Chernykh & Stepanov approach
% 1Aug18: Confirmed 4th order on forward iteration using soliton
    % hits machine limit at dz = 5e-4, nt = 512 (L=40)
% Quartic rate CONFIRMED for global iteration (4Aug18)

% domain: on real line truncated with periodic BCs at [-L,L), mapped
% to [-pi,pi)
% from z=0 to z=zeta, mapped to unit interval

% du/dz = i zeta(pi/L)^2 d^2u/dt^2 + i zeta |u|^2 u + zeta g
% u(0) = 0
% target at z=1: || sigma(u-v_zeta)||^2 = eta

% dg/dz = i zeta(pi/L)^2 d^2g/dt^2 + 2i zeta |u|^2 g - i zeta u^2 conj(g)
% g(1) = beta sigma^2 (u(1)-v_zeta)
