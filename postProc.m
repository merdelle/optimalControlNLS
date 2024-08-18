% script for processing optConNLS output
% run after running submitScript (or loading saved file)

betaVals = input.betaVals;
condVals = output.condVals;
gL2Vals = output.gL2Vals;
kappa = input.kappa;
zeta = input.zeta;
nt = input.nt; 
dz = input.dz; 
L = input.L;
betaKeeps = output.betaKeeps;
uKeeps = output.uKeeps;
gKeeps = output.gKeeps;
eta = input.eta;
sigma = input.sigma;

dt = 2*L/nt;
nz = round(1/dz);

ts = linspace(-L,L-dt,nt)';
zs = linspace(0,1,nz+1);

figure(1); plot(betaVals,condVals,betaVals,gL2Vals, ...
    betaVals,gL2Vals + kappa*condVals, ...
    betaVals,eta*ones(1,length(betaVals)))
xlabel('\beta')
legend({'$\|\sigma(u(\zeta)-v_{\zeta})\|^2$','$\|g\|^2$','J'}, ...
    'interpreter','latex');

for j = 1:length(betaKeeps),
    beta = betaKeeps(j);
    u = uKeeps(:,:,j);
    g = gKeeps(:,:,j);
    figure(2); mesh(zeta*zs,L/pi*ts,abs(u));
    xlabel('z'); ylabel('t'); zlabel('|u|');
    figure(3); mesh(zeta*zs,L/pi*ts,abs(g));
    xlabel('z'); ylabel('t'); zlabel('|g|');
    figure(4); plot(zeta*zs,atan2(imag(u(nt/2,:)),real(u(nt/2,:))),'o', ...
        zeta*zs,zeta*zs);
    xlabel('z');
    legend('arg(u)','A^2*z');
    figure(5); plot(zeta*zs,max(abs(u),[],1));
    xlabel('z'); ylabel('max(|u|)');
    figure(6); plot(L/pi*ts,abs(u(:,nz+1)).^2,L/pi*ts,sigma(ts));
    pause
end