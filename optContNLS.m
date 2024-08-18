function output = optContNLS(input);

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

% Parameters

L = input.L; %10;     % domain half-width
zeta = input.zeta; %10;  % propagation distance

kappa = input.kappa; %0;  % terminal cost factor
eta = input.eta; %0;    % terminal state margin
%dBeta = input.dBeta; %-0.01; % step size in beta to find zero of constraint
betaVals = input.betaVals; % array of betaVals for continuation
cCon = input.cCon; %.01; % convex sum factor to improve convergence

nt = input.nt; %1024; 
dz = input.dz; %0.01; 

uTol = input.uTol; %1e-5;    % tolerance in difference in u or diff in terminal cond
etaTol = input.etaTol; %1e-2;  % tolerance in constraint condition

sigma = input.sigma; %@(t) exp(-abs(t-tr)/tau); % window function
u0 = input.u0; % initial function for u

plotFlag = input.plotFlag; % set to 1 if want plots (not for parallel)
iterMax = input.iterMax;    % max num of iterations in u

% Initialization

dt = 2*L/nt;
nz = round(1/dz);

ts = linspace(-L,L-dt,nt)';
zs = linspace(0,1,nz+1);

g = zeros(nt,nz+1);   % initial costate
u = zeros(nt,nz+1); % initialize state
u(:,1) = u0(ts);

condVals = zeros(1,length(betaVals)); % values of terminal condition to plot
gL2Vals = condVals; % values of g norm to plot

betaKeeps = []; % keep track of beta values at special events
uKeeps = [];
gKeeps = [];

% Outer loop: adjust beta to achieve constraint

contFlag = 1; % turn off contFlag when done scanning through betas

eventFlag = 0; % flag to indicate when we reach target set

for j = 1:length(betaVals);
    
    beta = betaVals(j);
   
    % integrate u forward
    
 %   uOld = u;
    
    [u,g,termCond,errFlag] = optContNLS_singleBeta(input,beta,u,g);

    % useful output
    
    fprintf('beta = %6.4f, termCond = %6.4f\n',beta,termCond);
    gL2 = dz*dt*sum(trapz(abs(g).^2),2);
    
    condVals(j) = termCond;
    gL2Vals(j) = gL2;
    
    if beta == -kappa,
        fprintf('At beta=-kappa.\n');
        betaKeeps = [betaKeeps beta];
        uKeeps = cat(3,uKeeps,u);
        gKeeps = cat(3,gKeeps,g);
    end
    
    if errFlag == 1,
        fprintf('Failed to converge.\n');
        break
    elseif eventFlag == 1, % if already in eligible set
        if termCond > eta + etaTol,
            betaKeeps = [betaKeeps beta];
            uKeeps = cat(3,uKeeps,u);
            gKeeps = cat(3,gKeeps,g);
            fprintf('Left eligible set.  Terminating.\n');
            break
        end
    elseif termCond < eta + etaTol,
        betaKeeps = [betaKeeps beta];
        uKeeps = cat(3,uKeeps,u);
        gKeeps = cat(3,gKeeps,g);
        fprintf('Entering eligible set.\n');
        eventFlag = 1;
    end
    
%    iterNum

%     betaOld = beta;
%     if (iterNum == iterMax) || (max(max(abs(u)))>1e9) || isnan(termCond),
%         betaKeeps = [betaKeeps betaOld];
%         uKeeps = cat(3,uKeeps,uOld);
%         gKeeps = cat(3,gKeeps,gOld);
%         fprintf('Iteration failed.\n');
%         contFlag = 0;
%     elseif kappa == 0,
%         if termCond - eta < etaTol, % assumes termCond will start greater than eta
%             if abs(termCond - eta) < etaTol,
%                 betaKeeps = [betaKeeps beta];
%                 uKeeps = cat(3,uKeeps,u);
%                 gKeeps = cat(3,gKeeps,g);
%                 contFlag = 0;   % we are done, stop iterating
%             else % conduct bisection search in beta to refine target
%                 if contFlag == 2,
%                     if (termCond-eta)*(cvLeft-eta) > 0,
%                         betaLeft = beta;
%                         cvLeft = termCond;
%                     else
%                         betaRight = beta;
%                         cvRight = termCond;
%                     end
%                 else
%                     contFlag = 2; % signal entry into bisection method
%                     betaLeft = betaVals(length(betaVals)-1);
%                     cvLeft = condVals(length(betaVals)-1);
%                     betaRight = beta;
%                     cvRight = termCond;
%                 end
%                 beta = 1/2*(betaLeft+betaRight);
%                 g = beta/betaOld*g; % correct for beta to align termCond better
%             end
%         else
%             beta = beta + dBeta;    %
%         end
%     else % continue until beta = kappa or back on other side of termCond == 0
%         if (contFlag == 2) || (contFlag == 3),
%             if abs(beta + kappa) < etaTol,
%                 betaKeeps = [betaKeeps beta];
%                 uKeeps = cat(3,uKeeps,u);
%                 gKeeps = cat(3,gKeeps,g);
%                 contFlag = 0;   % done. found beta = -kappa
%             elseif (beta < -kappa) && (contFlag ~= 3),
%                 beta = -kappa;  % almost done. one more iteration
%             elseif (termCond - eta > etaTol)
%                 fprintf('Min must be at boundary.\n');
%                 contFlag = 0;   % Min must be at one of boundaries
%             else
%                 beta = beta + dBeta;
%             end
%         elseif termCond - eta < etaTol,
%             if beta < -kappa,
%                 fprintf('Warning: beta already past -kappa.\n');
%                 contFlag = 3; % Min must be at one of boundaries
%             else
%                 contFlag = 2;
%             end
%             betaKeeps = [betaKeeps beta];
%             uKeeps = cat(3,uKeeps,u);
%             gKeeps = cat(3,gKeeps,g);
%             beta = beta + dBeta;
%         else
%             beta = beta + dBeta;
%         end
%     end
end

if plotFlag == 1,
    figure(2); plot(betaVals,condVals,betaVals,gL2Vals, ...
        betaVals,gL2Vals + kappa*condVals, ...
        betaVals,eta*ones(1,length(betaVals)))
    xlabel('\beta')
    legend({'$\|\sigma(u(\zeta)-v_{\zeta})\|^2$','$\|g\|^2$','J'}, ...
        'interpreter','latex');
    figure(3); mesh(zeta*zs,ts,abs(u));
    xlabel('z'); ylabel('t'); zlabel('|u|');
    figure(4); mesh(zeta*zs,ts,abs(g));
    xlabel('z'); ylabel('t'); zlabel('|g|');
    figure(5); plot(zeta*zs,atan2(imag(u(nt/2,:)),real(u(nt/2,:))),'o', ...
        zeta*zs,atan2(imag(u(nt/2,:)),real(u(nt/2,:)))+2*pi,'o', ...
        zeta*zs,zeta*zs);
    figure(6); plot(zeta*zs,max(abs(u),[],1));
end

% define output

%output.betaVals = betaVals;
output.condVals = condVals;
output.gL2Vals = gL2Vals;
output.betaKeeps = betaKeeps;
output.uKeeps = uKeeps;
output.gKeeps = gKeeps;



