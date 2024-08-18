% test of whether iterative method for ODEs also exhibits 2nd order conv

nzs = 2.^10; %(1:10);
dzs = 1./nzs;
nN = length(nzs);
cCon = 1; %0.04;

mytol = 1e-13;

uFinals = zeros(1,nN);
gFirsts = zeros(1,nN);

errVals = zeros(1,nN-1);
    


for k = 1:nN,
    
    
    nz = nzs(k)
    dz = dzs(k);
    u = zeros(1,nz+1);
    u(1) = 1;
    g = zeros(1,nz+1);
    %zs = (0:dz:1);
    
    uDiff = 1e9;
    gDiff = 1e9;
    
    uDiffs = [];
    gDiffs = [];
    
    
%    while (uDiff+gDiff) > mytol,
 for jj = 1:2,
     
        uNew = u;
        for j = 1:nz,
            
            % advance  in u
            
            ul = uNew(j);
            gl = g(j);
            gr = g(j+1);
            if j == 1,
                grr = g(3);
                dgdz = 2i*abs(ul).^2.*gl -1i*ul.^2.*conj(gl);
                gll = -3/2*gl + 3*gr - 1/2*grr - 3*dz*dgdz;
            elseif j == nz,
                gll = g(nz-1);
                dgdz = 2i*abs(ul).^2.*gl -1i*ul.^2.*conj(gl);
                grr = 6*gr - 3*gl - 2*gll - 6*dz*dgdz;
            else
                gll = g(j-1);
                grr = g(j+2);
            end
            gm = 9/16*(gl+gr) - 1/16*(gll+grr);
            
            % RK4 step in u
            
            a = dz*(1i*abs(ul).^2.*ul + gl);
            u1 = ul + a/2;
            b = dz*(1i*abs(u1).^2.*u1 + gm);
            u2 = ul + b/2;
            c = dz*(1i*abs(u2).^2.*u2 + gm);
            u3 = ul + c;
            d = dz*(1i*abs(u3).^2.*u3 + gr);
            uNew(j+1) = ul + (a+2*b+2*c+d)/6;
            
        end
        
        uOld = u;
        u = cCon*uNew + (1-cCon)*u;
        
        uDiff = sqrt(dz*sum(abs(uNew-uOld).^2));
        
        % backstep in g

        gNew = g;
        gNew(nz+1) = 1-u(nz+1);
        
        for j = nz:-1:1,
            gr = gNew(j+1);
            ul = u(j);
            ur = u(j+1);
            if j == 1,
                urr = u(3);
                dudz = 1i*abs(ur).^2.*ur + gr;
                ull = 6*ul - 3*ur - 2*urr + 6*dz*dudz;
            elseif j == nz,
                ull = u(nz-1);
                dudz = 1i*abs(ur).^2.*ur + gr;
                urr = -3/2*ur + 3*ul - 1/2*ull + 3*dz*dudz;
            else
                ull = u(j-1);
                urr = u(j+2);
            end
            
            um = 9/16*(ul+ur) - 1/16*(ull+urr);
            
            % RK4 step in u
            
            a = -dz*(2i*abs(ur).^2.*gr - 1i*ur.^2.*conj(gr));
            g1 = gr + a/2;
            b = -dz*(2i*abs(um).^2.*g1 - 1i*um.^2.*conj(g1));
            g2 = gr + b/2;
            c = -dz*(2i*abs(um).^2.*g2 - 1i*um.^2.*conj(g2));
            g3 = gr + c;
            d = -dz*(2i*abs(ul).^2.*g3 - 1i*ul.^2.*conj(g3));
            gNew(j) = gr + (a+2*b+2*c+d)/6;
        end
        
        gDiff = sqrt(dz*sum(abs(gNew-g).^2));
        g = gNew;

%       [uDiff gDiff]
%        plot(zs,abs(u),zs,abs(g));
%        drawnow;

        if isnan(u(nz+1)+g(nz+1)),
            break
        end

        uDiffs = [uDiffs uDiff];
        gDiffs = [gDiffs gDiff];
        
    end

%     semilogy([uDiffs; gDiffs]','o');
%     pause
        

    uFinals(k) = u(nz+1);
    gFirsts(k) = g(1);
end

figure(1);
loglog(dzs(1:nN-1),abs(uFinals(1:nN-1)-uFinals(nN)),'o',...
    dzs(1:nN-1),abs(gFirsts(1:nN-1)-gFirsts(nN)),'s',...
    dzs(1:nN-1),dzs(1:nN-1),dzs(1:nN-1),dzs(1:nN-1).^4);
