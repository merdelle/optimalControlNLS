function [finalU,firstG,condVals] = getFinalUG(input,dz,ntVals,j);

% local function to address classification issues in parallelization
% For use only in convergence studies

nNt = length(ntVals);
ntMax = ntVals(nNt);
nt = ntVals(j);
nz = round(1/dz);

finalU = zeros(ntMax,1);
firstG = finalU;

tMesh = zeros(ntVals(nNt),1);
tMesh(1:2^(nNt-j):ntVals(nNt)-(2^(nNt-j)-1)) = 1;
tMesh = logical(tMesh);

input.nt = nt;
input.dz = dz;
output = optContNLS(input);
finalU(tMesh) = output.uKeeps(:,nz+1);
firstG(tMesh) = output.gKeeps(:,1);
condVals = output.condVals;

end


