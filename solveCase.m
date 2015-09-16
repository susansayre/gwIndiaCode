function [output] = solveCase(P,modelOpts,parameterSetID)

%% compute derived parameter values
P.h0 = P.landHeight - P.initialLift;

P.investInd = 1; P.gwInd = 2; P.wellInd = 1; P.levelInd = 2; P.yrsInd = 3;
%% solve the optimal management problem
model.func = 'optFunc';
model.discount = P.discount;
model.params = {P};

smin = [P.bottom 0];
smax = [P.landHeight P.maxWellCap];
n = [modelOpts.heightNodes modelOpts.capNodes];
fspace = fundefn('cheb',n,smin,smax);
snodes = funnode(fspace);
s = gridmake(snodes);
v = ones(size(s,1),1); x=[v v];
[c,s,v,x,resid] = dpsolve(model,fspace,s,v,x);
output.opt = v;
%extract and store necessary optimal management output

%% solve the adaptive expectations management problem
output.aeOut = aeSolve(P,modelOpts);

%% solve the "rational" expectations management problem
output.reOut = reSolve(P,modelOpts);

%% compare outputs and return key results

%% store details for later reference if needed
if ~exist('detailedOutput','dir')
    mkdir('detailedOutput')
end
if ~exist(['detailedOutput/' parameterSetID.timeStamp],'dir')
    mkdir(['detailedOutput/' parameterSetID.timeStamp])
end
save(['detailedOutput/' parameterSetID.timeStamp '/' parameterSetID.case ''])