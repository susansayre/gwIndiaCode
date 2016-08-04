function [output] = solveCase(P,modelOpts,parameterSetID)

%% compute derived parameter values
P.h0 = P.landHeight - P.initialLift;

P.investInd = 1; P.gwDugInd = 2; P.gwBoreInd = 3; P.sbInd = 1; P.levelInd = 2; P.yrsInd = 3;

P.inTruncProb = normcdf(P.maxInvestCost,P.investCostMean,P.investCostSD) - normcdf(P.minInvestCost,P.investCostMean,P.investCostSD);
P.probBelow = normcdf(P.minInvestCost,P.investCostMean,P.investCostSD);
P.lowCost = norminv(P.shrBore0*P.inTruncProb + P.probBelow,P.investCostMean,P.investCostSD);

%% solve the optimal management problem
model.func = 'optFunc';
model.discount = P.discount;
%model.horizon = 25;

smin(P.levelInd) = P.bottom;
%smin(P.levelInd) = 10;
smax(P.levelInd) = P.landHeight;
smin(P.sbInd) = P.shrBore0;
smax(P.sbInd) = 1;
n = [modelOpts.capNodes modelOpts.heightNodes];
model.ds = numel(n);
fspace = fundefn('spli',n,smin,smax);
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
costBore = P.electricityBore*lifts;
costDug = P.electricityDug*lifts;
P.dx = 3;

states = length(s);
[lb,ub] = feval(model.func,'b',s,ones(states,3),[],[],[],[],P);

xGuess(:,P.gwDugInd) = max((P.idDugInt - costDug)./P.idDugSlope,0);
xGuess(:,P.gwBoreInd) = max((P.idBoreInt - costBore)./P.idBoreSlope,0);

nb = netBen(xGuess(:,P.gwDugInd),xGuess(:,P.gwBoreInd),xGuess(:,P.investInd),s(:,P.levelInd),s(:,P.sbInd),P);
deltaBen = nb.bore - nb.dug;
newShr = max(s(:,P.sbInd),normcdf(deltaBen*P.discount/(1-P.discount),P.investCostMean,P.investCostSD));
xGuess(:,P.investInd) = newShr-s(:,P.sbInd);

xGuess = max(lb,min(xGuess,ub));
if isfield(model,'horizon'); finite = 1; else, finite = 0; end

vGuess = feval(model.func,'f',s,xGuess,[],[],[],[],P)/(1-P.discount);
if ~finite
    %xGuess = .9*xGuess;
    %xGuess(:,P.investInd) = min(.5*ub(:,P.investInd),newShr-s(:,P.sbInd));
end    

gLB = min(lb);
gUB = max(ub);

model.dx = size(xGuess,2);
P.dx = model.dx;
model.params = {P};

optset('dpsolve','algorithm',modelOpts.algorithm);
optset('dpsolve','maxit',modelOpts.maxit);
optset('dpsolve','showiters',1);
optset('dpsolve','tol',modelOpts.vtol);
optset('dpsolve','nres',1);
optset('dpsolve','maxitncp',200);
optset('dpsolve_vmax','maxbacksteps',0);
optset('dpsolve_vmax','maxit',50);
optset('dpsolve_vmax','lcpmethod','minmax');

save beforeDPsolve

%figure()
%begin by solving on finite grid of x values to get close
modelApx = model;
xApproxNodes = [20 20 20];
for ii=1:numel(xApproxNodes)
    stepSize = (gUB(ii) - gLB(ii))/(xApproxNodes(ii)-1);
    xNodes{ii} = (gLB(ii):stepSize:gUB(ii))';
end
modelApx.X = gridmake(xNodes);  
optset('dpsolve','tol',.1);
[~,~,v,x] = dpsolve3(modelApx,fspace,vGuess,xGuess);

optset('dpsolve','tol',modelOpts.vtol);
[c,sr,vr,xr,resid] = dpsolve3(model,fspace,v,x);
max(abs(resid./vr))
%[c,scoord,v,x] = dpsolve(model,fspace,s,vGuess,xGuess);

%output.opt.converged = exf;

%[shares,levels] = ndgrid(scoord{P.sbInd},scoord{P.levelInd});
s0(P.sbInd) = P.shrBore0;
s0(P.levelInd) = P.h0;
[ssim,xsim] = dpsimul3(model,fspace,modelOpts.minT,s0,1,sr,vr,xr);

%extract and store necessary optimal management output
output.opt.vFunc = vr;
output.opt.val = funeval(c,fspace,s0);
output.opt.statePath = squeeze(ssim);
output.opt.controlPath = squeeze(xsim);
optNb = netBen(output.opt.controlPath(:,P.gwDugInd),output.opt.controlPath(:,P.gwBoreInd),output.opt.controlPath(:,P.investInd),output.opt.statePath(:,P.levelInd),output.opt.statePath(:,P.sbInd),P);
output.opt.valPath = [optNb.dug optNb.bore optNb.all];
output.opt.optVal = (P.discount.^(0:length(output.opt.valPath)-1))*output.opt.valPath(:,3);

% save beforeAe
% % solve the adaptive expectations management problem
% output.aeOut = aeSolve(P,modelOpts);
% output.aeOut.pgain = (output.opt.val - output.aeOut.aeVal)/output.aeOut.aeVal;
% 
% if ~isfield(modelOpts,'figureVisible')
%     modelOpts.figureVisible = 'off';
% end
% 
%figure('Visible',modelOpts.figureVisible)
xTitles = {'Investment','(Traditional)','(Modern)'};
sTitles = {'Share in Modern Agriculture','Pumping Lift'};
sYlabel = {'%','meters'};
% 
 save beforeRe
% % solve the "rational" expectations management problem
%  output.reOut = reSolve(P,modelOpts,output.aeOut.statePath(:,2));
%  output.reOut.pgain = (output.opt.val - output.reOut.reVal)/output.reOut.reVal;

 output.reOut = reSolve(P,modelOpts,output.opt.statePath(:,2));
 %output.reOut = reSolve(P,modelOpts,P.h0*ones(100,1));
 output.reOut.pgain = (output.opt.val - output.reOut.reVal)/output.reOut.reVal;

 yrs = min([length(output.reOut.controlPath) length(output.opt.controlPath)]);

 byFarmFig = figure();
 for ii=2:3; 
    subplot(2,2,ii-1); 
    hold on; 
    plot(output.reOut.controlPath(1:yrs,ii),'-.'); 
    plot(output.opt.controlPath(1:yrs,ii)); 
    title(['Water ' xTitles{ii}]);
    subplot(2,2,ii+1);
    hold on;
    plot(output.reOut.valPath(1:yrs,ii-1),'-.'); 
    plot(output.opt.valPath(1:yrs,ii-1));
    title(['Current Payoff ' xTitles{ii}])
   
end; 

totFig = figure();
for ii=1:2; 
    subplot(2,2,ii); 
    if ii==2
        constant = P.landHeight; slope = -1;
    else
        constant = 0; slope=1;
    end
    hold on; 
    plot(constant+slope*output.reOut.statePath(1:yrs,ii),'-.'); 
    plot(constant+slope*output.opt.statePath(1:yrs,ii)); 
    title(sTitles{ii});
    ylabel(sYlabel{ii});
end;

subplot(2,2,3); 
hold on; 
plot(output.reOut.statePath(1:yrs,1).*output.reOut.controlPath(1:yrs,3)+(1-output.reOut.statePath(1:yrs,1)).*output.reOut.controlPath(1:yrs,2),'-.'); 
plot(output.opt.statePath(1:yrs,1).*output.opt.controlPath(1:yrs,3)+(1-output.opt.statePath(1:yrs,1)).*output.opt.controlPath(1:yrs,2)); 
title('Total Water Use');

subplot(2,2,4); 
hold on; 
plot(output.reOut.valPath(1:yrs,3),'-.'); 
plot(output.opt.valPath(1:yrs,3)); 
title('Total Current Value');

%% compare outputs and return key results

%% store details for later reference if needed
if ~exist('detailedOutput','dir')
    mkdir('detailedOutput')
end
if ~exist(['detailedOutput/' parameterSetID{1}],'dir')
    mkdir(['detailedOutput/' parameterSetID{1}])
end
saveas(totFig,fullfile('detailedOutput',parameterSetID{1},[parameterSetID{2} '_AggPaths']),'epsc')
saveas(byFarmFig,fullfile('detailedOutput',parameterSetID{1},[parameterSetID{2} '_FarmPaths']),'epsc')
clear totFig byFarmFig
save(['detailedOutput/' parameterSetID{1} '/' parameterSetID{2} ''])