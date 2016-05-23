function [output] = solveCase(P,modelOpts,parameterSetID)

%% compute derived parameter values
P.h0 = P.landHeight - P.initialLift;

P.investInd = 1; P.gwDugInd = 2; P.gwBoreInd = 3; P.sbInd = 1; P.levelInd = 2; P.yrsInd = 3;

%Set inverse demand curve parameters for traditional farm
P.idDugInt = 2*P.dDugInt; %normalizes values in terms of max benefit of free water on traditional farm
P.idDugSlope = P.idDugInt/P.dDugInt;

P.idBoreInt = 2*P.boreVInc/(P.dDugInt*P.boreQInc); %determines vertical intercept that guarantees max value = boreVInc * max value on trad farm
P.idBoreSlope = P.idBoreInt/(P.dDugInt*P.boreQInc); %determines the slope consistent with pInt = P.idBoreInt and qInt = P.dDugInt*P.boreQInc

%per unit cost function for dug wells has form C(h) = ae^(bh). Set so its
%value is equal to the demand intercept at maxDepthDug and its slope is
%equal to slopeMaxDepth at this same value.
P.costDug_a = P.dDugInt*exp(-P.maxDepthDug*P.slopeMaxDepth/P.dDugInt); P.costDug_b = P.slopeMaxDepth/P.dDugInt;

%add a cost as we near the bottom of the aquifer that ensures it is
%uneconomic to pump water at the bottom of the aquifer.
P.costBore_a = log(P.idBoreInt);
P.costBore_b = P.slopeBoreZero*exp(-P.costBore_a);
P.costBore_a = 0;

P.inTruncProb = normcdf(P.maxInvestCost,P.investCostMean,P.investCostSD) - normcdf(P.minInvestCost,P.investCostMean,P.investCostSD);
P.probBelow = normcdf(P.minInvestCost,P.investCostMean,P.investCostSD);
P.lowCost = norminv(P.shrBore0*P.inTruncProb + P.probBelow,P.investCostMean,P.investCostSD);

%% solve the optimal management problem
model.func = 'optFunc';
model.discount = P.discount;
model.params = {P};
%model.horizon = 25;

smin(P.levelInd) = P.bottom;
%smin(P.levelInd) = 10;
smax(P.levelInd) = P.landHeight;
smin(P.sbInd) = P.shrBore0;
smax(P.sbInd) = 1;

n = [modelOpts.capNodes modelOpts.heightNodes];
fspace = fundefn('lin',n,smin,smax);
minHeightDug = P.landHeight - P.maxDepthDug;
%Add a dense grid of points close to maxDepthDug
% heightBreaks = fspace.parms{P.levelInd}{1};
% fspace.parms{P.levelInd}{1} = unique(sort([heightBreaks; (minHeightDug-1:.1:minHeightDug+1)']));
% fspace.parms{P.levelInd}{2} = 0;
% fspace.n(P.levelInd) = length(fspace.parms{P.levelInd}{1});
% 
% %Add a dense grid of points near begining and end of investment intervals
% capBreaks = fspace.parms{P.sbInd}{1};
% fspace.parms{P.sbInd}{1} = unique(sort([capBreaks; (P.shrBore0+.01:.01:P.shrBore0+.2)'; (.95:.01:.99)']));
% fspace.parms{P.sbInd}{2} = 0;
% fspace.n(P.sbInd) = length(fspace.parms{P.sbInd}{1});

snodes = funnode(fspace);
s = gridmake(snodes);

%make guess at the optimal actions
lifts = P.landHeight - s(:,P.levelInd);
costBore = P.electricity*lifts;
costDug = P.costDug_a*exp(P.costDug_b*lifts);

states = length(s);
[lb,ub] = feval(model.func,'b',s,ones(states,3),[],P);

xGuess(:,P.gwDugInd) = max((P.idDugInt - costDug)./P.idDugSlope,0);
xGuess(:,P.gwBoreInd) = max((P.idBoreInt - costBore)./P.idBoreSlope,0);

nb = netBen(xGuess(:,P.gwDugInd),xGuess(:,P.gwBoreInd),xGuess(:,P.investInd),s(:,P.levelInd),s(:,P.sbInd),P);
deltaBen = nb.bore - nb.dug;
newShr = max(s(:,P.sbInd),normcdf(deltaBen*P.discount/(1-P.discount),P.investCostMean,P.investCostSD));

if isfield(model,'horizon'); finite = 1; else, finite = 0; end

vGuess = feval(model.func,'f',s,xGuess,[],P)/(1-P.discount);
if ~finite
    %xGuess = .9*xGuess;
    %xGuess(:,P.investInd) = min(.5*ub(:,P.investInd),newShr-s(:,P.sbInd));
end    

optset('dpsolve','algorithm',modelOpts.algorithm);
optset('dpsolve','maxit',modelOpts.maxit);
optset('dpsolve','showiters',1);
optset('dpsolve','tol',modelOpts.vtol);
optset('dpsolve','nres',1);
optset('dpsolve_vmax','maxbacksteps',0);
optset('dpsolve_vmax','maxit',50);
optset('dpsolve_vmax','lcpmethod','minmax');

%figure()
[c,scoord,v,x,resid,exf] = dpsolveS2(model,fspace,s,vGuess,xGuess);
%[c,scoord,v,x] = dpsolve(model,fspace,s,vGuess,xGuess);

output.opt.converged = exf;

[shares,levels] = ndgrid(scoord{P.sbInd},scoord{P.levelInd});
s0(P.sbInd) = P.shrBore0;
s0(P.levelInd) = P.h0;
[ssim,xsim] = dpsimul(model,s0,modelOpts.minT,scoord,x);
close

output.opt.vFunc = v;
output.opt.val = funeval(c,fspace,s0);
output.opt.statePath = squeeze(ssim)';
output.opt.controlPath = squeeze(xsim)';
output.opt.valPath = feval(model.func,'f',output.opt.statePath,output.opt.controlPath,[],P);
output.opt.optVal = (P.discount.^(0:length(output.opt.valPath)-1))*output.opt.valPath;
%extract and store necessary optimal management output

% solve the adaptive expectations management problem
output.aeOut = aeSolve(P,modelOpts);
output.aeOut.pgain = (output.opt.val - output.aeOut.aeVal)/output.aeOut.aeVal;

if ~isfield(modelOpts,'figureVisible')
    modelOpts.figureVisible = 'off';
end

figure('Visible',modelOpts.figureVisible)
xTitles = {'Investment','Water (Traditional)','Water (Modern)'};
sTitles = {'Share in Modern Agriculture','Pumping Lift'};
sYlabel = {'%','meters'};

% solve the "rational" expectations management problem
 output.reOut = reSolve(P,modelOpts);
 output.reOut.pgain = (output.opt.val - output.reOut.reVal)/output.reOut.reVal;
 
 yrs = min([length(output.reOut.controlPath) length(output.aeOut.controlPath) length(output.opt.controlPath)]);

 for ii=2:3; 
    subplot(2,2,ii+1); 
    plot(output.aeOut.controlPath(1:yrs,ii),'--'); 
    hold on; 
    plot(output.reOut.controlPath(1:yrs,ii),'-.'); 
    plot(output.opt.controlPath(1:yrs,ii)); 
    title(xTitles{ii});
    legend('Common Property AE','Common Property RE','Optimal Management')
end; 

for ii=1:2; 
    subplot(2,2,ii); 
    if ii==2
        constant = P.landHeight; slope = -1;
    else
        constant = 0; slope=1;
    end
    plot(constant+slope*output.aeOut.statePath(1:yrs,ii),'--'); 
    hold on; 
    plot(constant+slope*output.reOut.statePath(1:yrs,ii),'-.'); 
    plot(constant+slope*output.opt.statePath(1:yrs,ii)); 
    title(sTitles{ii});
    ylabel(sYlabel{ii});
    legend('Common Property AE','Common Property RE','Optimal Management')
end;



%% compare outputs and return key results

%% store details for later reference if needed
if ~exist('detailedOutput','dir')
    mkdir('detailedOutput')
end
if ~exist(['detailedOutput/' parameterSetID.timeStamp],'dir')
    mkdir(['detailedOutput/' parameterSetID.timeStamp])
end
save(['detailedOutput/' parameterSetID.timeStamp '/' parameterSetID.case ''])
saveas(gcf,fullfile('detailedOutput',parameterSetID.timeStamp,[parameterSetID.case '_paths']),'epsc')
