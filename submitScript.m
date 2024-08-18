%clear all; close all;

% submit script for function optimalControlNLS

input.L = 40;     % domain half-width, should be at least 40 for 1-soliton to avoid boundary effects
input.zeta = 5;  % propagation distance

input.kappa = 1;  % terminal cost factor
input.eta = 0.1;    % terminal state margin
%input.dBeta = -0.01; % step size in beta to find zero of constraint
input.cCon = 0.05; % convex sum factor to improve convergence
input.betaVals = (0:-0.5:-10);

input.nt = 512; %2048; %1024;
input.dz = 1e-4; %5e-4; %0.01;

input.uTol = 1e-6;    % tolerance in difference in u or diff in terminal cond
input.etaTol = 1e-3;  % tolerance in constraint condition

tr = 1;     % window centre
tau = 4;    % window width
input.sigma = @(t) exp(-(t-tr).^2/tau^2); % window function
%input.sigma = @(t) exp(-abs(t-tr)/tau); % window function
%input.sigma = @(t) ones(size(t));

A = 0;
phi = 0;
tFinal = 0; %2;
Om = 0;
input.vTarg = @(t) sqrt(2)*A*sech(A*(t-tFinal)).*exp(1i*A^2*input.zeta+1i*phi+1i*Om*t); % target function for u
%input.vTarg = @(t) ones(size(t));

A0 = 1;
%input.u0 = sqrt(2)*A0*sech(A0*ts);   % initial condition for u
input.u0 = @(t) sqrt(2)*A0*sech(A0*t);
%input.u0 = @(t) ones(size(t));

input.plotFlag = 0; % set to 1 if want plots (not for parallel)
input.iterMax = 1000;

%% Describe runs

runDescript = 'Killing soliton.  Medium shifted window.\n';
output = optContNLS(input);
save

%% Convergence studies using last slice of u, and first slices of g

% runDescript = 'Convergence study in dz for full iteration growing one-soliton.\n';
% 
% ntVals = 2.^(4:13);
% dz = input.dz;
% nz = round(1/dz);
% %dzVals = 1./(2.^(4:12)); %(6:10)); %15));
% nt = input.nt;
% L = input.L;
% dt = 2*L/nt;
% nN = length(dzVals); %ntVals);
% errVals = zeros(1,nN-1);
% finalUs = zeros(nt,nN); %ntVals(nN),nN);
% firstGs = finalUs;
% allCondVals = zeros(length(input.betaVals),nN);
% 
% kv = [0:nt/2-1 -nt/2:-1]';
% 
% ts = linspace(-L,L-dt,nt)';
% vTarg = input.vTarg(ts);
% beta = input.betaVals;
% zeta = input.zeta;
% 
% %for k = 1:nN,
%  
% parfor labindex = 1:nN,
%     dz = dzVals(labindex);
% 
% %    nz = round(1/dz);
% %     for j = 1:nz+1,
% %         u(:,j) = sqrt(2)*sech(ts).*exp(1i*zeta*dz*(j-1));
% %     end
% %     g(:,nz+1) = u(:,nz+1);
% 
% %     for j = nz:-1:1,
% %         uRight = sqrt(zeta)*u(:,j+1);
% %         uLeft = sqrt(zeta)*u(:,j);
% %         if j == nz,
% %             % build O(dz^2) approx for u(nz+2) using du/dz
% %            dudz = zeta*(ifft(-1i*(pi/L)^2*kv.^2.*fft(u(:,nz+1))) + ...
% %                1i*abs(u(:,nz+1)).^2.*u(:,nz+1)); % + g(:,nz+1));
% %            uRR = -3/2*uRight + 3*dz*sqrt(zeta)*dudz + 3*sqrt(zeta)*u(:,nz) - ...
% %                1/2*sqrt(zeta)*u(:,nz-1);
% %         else
% %             uRR = sqrt(zeta)*u(:,j+2);
% %         end
% %         if j == 1,
% %             % build O(dz^2) approx for u(0) using du/dz
% %            dudz = zeta*(ifft(-1i*(pi/L)^2*kv.^2.*fft(u(:,1))) + ...
% %                1i*abs(u(:,1)).^2.*u(:,1)); % + g(:,1));
% %            uLL = -3/2*uLeft - 3*dz*sqrt(zeta)*dudz + 3*sqrt(zeta)*u(:,2) - ...
% %                1/2*sqrt(zeta)*u(:,3);
% %         else
% %             uLL = sqrt(zeta)*u(:,j-1);
% %         end
% %         gRight = g(:,j+1);
% %         gLeft = lnls_RK4Istep(gRight,uLL,uLeft,uRight,uRR,zeta*(pi/L)^2,nt,dz);
% %         g(:,j) = gLeft;
% %     end
% 
% %    errVals(k) = sqrt(dt*trapz(abs(uNew(:,nz+1)-vTarg).^2));
% 
%     
%     [finalU,firstG,condVals] = getFinalUG(input,dz,nt,1); %ntVals,j);
% 
%     finalUs(:,labindex) = finalU;
%     firstGs(:,labindex) = firstG;
%     allCondVals(:,labindex) = condVals;
% end
% 
% %loglog(dzVals,errVals,'o',dzVals,dzVals.^4)
% 
% for k = 1:nN-1,
%     
% %     tMesh = zeros(nN,1);
% %     tMesh(1:2^(nN-k):2^nN-(2^(nN-k)-1)) = 1;
% %     tMesh = logical(tMesh);
% 
%     tMesh = logical(ones(nt,1));
% 
% %    dt = 2*input.L/ntVals(k);
%    errVals(k) = sqrt(dt*(sum(abs(finalUs(tMesh,k)-finalUs(tMesh,nN)).^2))+ ...
%        sum(abs(firstGs(tMesh,k)-firstGs(tMesh,nN)).^2));
% 
% %    errVals(k) = sqrt(dt*trapz(abs(firstGs(:,k)-firstGs(:,nN)).^2));
% 
% end
% 
% if input.plotFlag == 1,
%     loglog(dzVals(1:nN-1),errVals,'o',dzVals,10^3*dzVals.^4)
%     xlabel('dz'); ylabel('RMS error');
%     %title('Convergence in dz for back-evolution of costate');
%     title('Convergence in dz using growth of 1-soliton');
% end
% 
