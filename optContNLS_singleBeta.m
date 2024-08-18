function [u,g,termCond,errFlag] = optContNLS_singleBeta(input,beta,u,g);

%
% component of optContNLS that performs iteration for fixed beta
%
% NOTE: to minimize number of parameters sent to RK4 functions, the loops
% below that call these functions to integrate forward and backward use the
% following reparameterizations:
%   forward (nls_RK4Istep): z -> z/zeta, u -> sqrt(zeta)u, g ->
%   (zeta)^(3/2) g
%
%   backward (lnls_RK4Istep): z -> z/zeta, u -> sqrt(zeta)u, g -> g
%

L = input.L; %10;     % domain half-width
zeta = input.zeta; %10;  % propagation distance

cCon = input.cCon; %.01; % convex sum factor to improve convergence

nt = input.nt; %1024; 
dz = input.dz; %0.01; 

uTol = input.uTol; %1e-5;    % tolerance in difference in u or diff in terminal cond

sigma = input.sigma; %@(t) exp(-abs(t-tr)/tau); % window function

vTarg = input.vTarg; %sqrt(2)*A*sech(A*ts)*exp(1i*A^2*zeta+1i);   % target function for u

plotFlag = input.plotFlag; % set to 1 if want plots (not for parallel)
iterMax = input.iterMax;    % max num of iterations in u

% Initialization

dt = 2*L/nt;
nz = round(1/dz);

ts = linspace(-L,L-dt,nt)';
zs = linspace(0,1,nz+1);
kv = [0:nt/2-1 -nt/2:-1]';

iterNum = 0;

termCondOld = 1e9;
termCond = 0;
uDiff = 1e9;

%while (uDiff > uTol) && (iterNum < iterMax),
while (abs(termCondOld-termCond) > uTol) && (iterNum < iterMax),
    
    uNew = u;
    for j = 1:nz,
        uLeft = sqrt(zeta)*uNew(:,j);
        gLeft = zeta^(3/2)*g(:,j);
        gRight = zeta^(3/2)*g(:,j+1);
        %
        if j == 1,
            gRR = zeta^(3/2)*g(:,3);
            % build O(dz^2) approx for g(-1) using du/dz
            dgdz = zeta*(ifft(-1i*(pi/L)^2*kv.^2.*fft(g(:,1))) + ...
                2i*abs(uNew(:,1)).^2.*g(:,1) - 1i*uNew(:,1).^2.*conj(g(:,1)));
            gLL = -3/2*gLeft - 3*zeta^(3/2)*dz*dgdz + 3*gRight - 1/2*gRR;
        elseif j == nz,
            gLL = zeta^(3/2)*g(:,nz-1);
            % build O(dz^2) approx for g(nz+2) using du/dz
            % Make sure not to use u(:,nz+1) in construction!
            dgdz = zeta*(ifft(-1i*(pi/L)^2*kv.^2.*fft(g(:,nz))) + ...
                2i*abs(uNew(:,nz)).^2.*g(:,nz) - 1i*uNew(:,nz).^2.*conj(g(:,nz)));
            gRR = 6*gRight - 3*gLeft - 2*gLL - 6*zeta^(3/2)*dz*dgdz;
        else
            gLL = zeta^(3/2)*g(:,j-1);
            gRR = zeta^(3/2)*g(:,j+2);
        end
               
        uRight = nls_RK4Istep(uLeft,gLL,gLeft,gRight,gRR,zeta*(pi/L)^2,nt,dz);
        uNew(:,j+1) = uRight/sqrt(zeta);
    end
    
    uOld = u;
    u = cCon*uNew + (1-cCon)*u; % to improve convergence
    
    uDiff = sqrt(dz*dt*sum(trapz(abs(uNew-uOld).^2,2)));
    
    gNew(:,nz+1) = beta*sigma(ts).^2.*(u(:,nz+1)-vTarg(ts));
    
    for j = nz:-1:1,
        uRight = sqrt(zeta)*u(:,j+1);
        uLeft = sqrt(zeta)*u(:,j);
        %
        if j == nz,
            uLL = sqrt(zeta)*u(:,nz-1);
            % build O(dz^2) approx for u(nz+2) using du/dz
            dudz = zeta*(ifft(-1i*(pi/L)^2*kv.^2.*fft(u(:,nz+1))) + ...
                1i*abs(u(:,nz+1)).^2.*u(:,nz+1) + gNew(:,nz+1));
            uRR = -3/2*uRight + 3*dz*sqrt(zeta)*dudz + 3*uLeft - 1/2*uLL;
        elseif j == 1,
            uRR = sqrt(zeta)*u(:,3);
            % build O(dz^2) approx for u(0) using du/dz
            % Make sure not to use g(:,1) in construction!
            dudz = zeta*(ifft(-1i*(pi/L)^2*kv.^2.*fft(u(:,2))) + ...
                1i*abs(u(:,2)).^2.*u(:,2) + gNew(:,2));
            uLL = 6*uLeft - 3*uRight - 2*uRR + 6*sqrt(zeta)*dz*dudz;
        else
            uLL = sqrt(zeta)*u(:,j-1);
            uRR = sqrt(zeta)*u(:,j+2);
        end
%
        gRight = gNew(:,j+1);
        gLeft = lnls_RK4Istep(gRight,uLL,uLeft,uRight,uRR,zeta*(pi/L)^2,nt,dz);
        gNew(:,j) = gLeft;
    end
    
    gDiff = sqrt(dz*dt*sum(trapz(abs(gNew-g).^2,2)));
    
    g = gNew;
    
    if plotFlag == 1,
        figure(1);
        h = plot(ts,abs(u(:,nz+1)),ts,abs(vTarg(ts)),ts,abs(g(:,nz+1)));
        axis([-L L -sqrt(2) sqrt(2)]);
        %            h(1).YData = real(u(:,nz+1));
        %            h(3).YData = real(g(:,nz+1));
        drawnow;
    end
    termCondOld = termCond;
    termCond = dt*trapz(abs(sigma(ts).*(u(:,nz+1)-vTarg(ts))).^2);
%    termCond-termCondOld
    
    iterNum = iterNum + 1;
end

errFlag = 0;
if (iterNum == iterMax) || (max(max(abs(u)))>1e9) || (isnan(termCond)),
    errFlag = 1;
end

