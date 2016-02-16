function [output] = solveCase(P,modelOpts,parameterSetID)

%% compute derived parameter values
P.h0 = P.landHeight - P.initialLift;

P.investInd = 1; P.gwInd = 2; P.wellInd = 1; P.levelInd = 2; P.yrsInd = 3;
%% solve the optimal management problem
model.func = 'optFunc';
model.discount = P.discount;
model.params = {P};

smin(P.levelInd) = P.bottom;
smax(P.levelInd) = P.landHeight;
smin(P.wellInd) = 0;
smax(P.wellInd) = P.maxWellCap;

n = [modelOpts.heightNodes modelOpts.capNodes];
fspace = fundefn('spli',n,smin,smax);

%force a breakpoint at max depth for dug well and add a repeated breakpoint
%at 0 well capital
minLevel = P.landHeight - P.maxDepthNoCap;
tempBreaks = fspace.parms{P.levelInd}{1};
tempBreaks = sort([tempBreaks; minLevel; minLevel-.1]);
fspace.parms{P.levelInd}{1} = tempBreaks;
fspace.parms{P.wellInd}{1} = sort([.001; fspace.parms{P.wellInd}{1}]);
snodes = funnode(fspace);
s = gridmake(snodes);

%make guess at the optimal actions
extraZeroRow = max(find(snodes{1}==0));
extraLandHeightCol  = min(find(snodes{2}==P.landHeight));
rows = numel(snodes{1}); cols = numel(snodes{2});
sHat = s;
sHat(sub2ind([rows cols 2],extraZeroRow*ones(1,cols),1:cols),ones(1,cols))=eps;
sHat(sub2ind([rows cols 2],1:rows,extraLandHeightCol*ones(1,rows),2*ones(1,rows))) = P.landHeight-eps;
lifts = P.landHeight - sHat(:,P.levelInd);
costs = P.electricity*sHat(:,P.wellInd).*lifts;

[lb,ub] = feval(model.func,'b',sHat,ones(size(s)),[],P);
xGuess(:,P.gwInd) = min(max(lb(:,P.gwInd),P.gwdIntercept - P.gwdSlope*costs),ub(:,P.gwInd));
xGuess(:,P.investInd) = .9*ub(:,P.investInd);

xGuess = .9*xGuess;

vGuess = feval(model.func,'f',sHat,xGuess,[],P);
optset('dpsolve','algorithm','newton');

[c,scoord,v,x,resid] = dpsolve(model,fspace,s,vGuess,xGuess);
s0(P.wellInd) = 0;
s0(P.levelInd) = P.h0;
[ssim,xsim] = dpsimul(model,s0,modelOpts.minT,scoord,xOut);
keyboard
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