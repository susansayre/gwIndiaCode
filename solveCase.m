function [output] = solveCase(P,modelOpts,parameterSetID)

%% compute derived parameter values
P.h0 = P.landHeight - P.initialLift;

P.investInd = 1; P.gwDugInd = 2; P.gwBoreInd = 3; P.sbInd = 1; P.levelInd = 2; P.yrsInd = 3;

%per unit cost function for dug wells has form C(h) = ae^(bh). Set so its
%value is equal to the demand intercept at maxDepthDug and its slope is
%equal to slopeMaxDepth at this same value.
P.costDug_a = P.dDugInt*exp(-P.maxDepthDug*P.slopeMaxDepth/P.dDugInt); P.costDug_b = log(P.dDugInt/P.costDug_a)/P.maxDepthDug;

%scaling parameters
P.valScale = 1/2*max(P.dDugInt^2/P.dDugSlope,P.dBoreInt^2/P.dBoreSlope);
%% solve the optimal management problem
model.func = 'optFunc';
model.discount = P.discount;
model.params = {P};

smin(P.levelInd) = P.bottom;
smax(P.levelInd) = P.landHeight;
smin(P.sbInd) = 0.01;
smax(P.sbInd) = 1;

n = [modelOpts.heightNodes modelOpts.capNodes];
fspace = fundefn('spli',n,smin,smax);

snodes = funnode(fspace);
s = gridmake(snodes);

%make guess at the optimal actions
lifts = P.landHeight - s(:,P.levelInd);
costBore = P.electricity*lifts;
costDug = P.costDug_a*exp(P.costDug_b*lifts);

states = length(s);
[lb,ub] = feval(model.func,'b',s,ones(states,3),[],P);
xGuess(:,P.gwDugInd) = max((P.dDugInt - costDug)./P.dDugSlope,0);
xGuess(:,P.gwBoreInd) = max((P.dBoreInt - costBore)./P.dBoreSlope,0);
nb = netBen(xGuess(:,P.gwDugInd),xGuess(:,P.gwBoreInd),xGuess(:,P.investInd),s(:,P.levelInd),s(:,P.sbInd),P);
[nb.dug nb.bore]
xGuess(:,P.investInd) = .5*ub(:,P.investInd);

xGuess = .9*xGuess;

vGuess = feval(model.func,'f',s,xGuess,[],P);
optset('dpsolve','algorithm','newton');
optset('dpsolve','maxit',2000);

[c,scoord,v,x,resid] = dpsolve(model,fspace,s,vGuess,xGuess);
[levels,shares] = ndgrid(scoord{P.levelInd},scoord{P.sbInd});
s0(P.sbInd) = 0;
s0(P.levelInd) = P.h0;
[ssim,xsim] = dpsimul(model,s0,modelOpts.minT,scoord,x);
figure()
subplot(2,1,1); plot(squeeze(ssim(1,1,:)));subplot(2,1,2); plot(squeeze(ssim(1,2,:)));

figure()
for ii=1:3; subplot(3,1,ii); plot(squeeze(xsim(1,ii,:))); end;

output.opt = v;
%extract and store necessary optimal management output

% solve the adaptive expectations management problem
output.aeOut = aeSolve(P,modelOpts);

figure()
for ii=1:3; subplot(3,1,ii); plot(output.aeOut.controlPath(:,ii)); end;
figure()
for ii=1:2; subplot(2,1,ii); plot(output.aeOut.statePath(:,ii)); end;

% solve the "rational" expectations management problem
 save beforeRe
 return
% output.reOut = reSolve(P,modelOpts,output.aeOut.statePath(:,P.levelInd),output.aeOut.controlPath);

%% compare outputs and return key results

%% store details for later reference if needed
if ~exist('detailedOutput','dir')
    mkdir('detailedOutput')
end
if ~exist(['detailedOutput/' parameterSetID.timeStamp],'dir')
    mkdir(['detailedOutput/' parameterSetID.timeStamp])
end
save(['detailedOutput/' parameterSetID.timeStamp '/' parameterSetID.case ''])