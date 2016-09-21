function [output] = solveCase(P,modelOpts,parameterSetID)

%check whether the case is investment limited
if strfind(parameterSetID{2},'niP0')
    solveType = 'niP0';
elseif strfind(parameterSetID{2},'niPEnd')
    solveType = 'niPEnd';
elseif strfind(parameterSetID{2},'niPAvg')
    solveType = 'niPAvg';
else
    solveType = 'regular';
end

if strcmp(solveType,'regular')
    %% compute derived parameter values
    P.h0 = P.landHeight - P.initialLift;

    P.investInd = 1; P.gwDugInd = 2; P.gwBoreInd = 3; P.sbInd = 1; P.levelInd = 2; P.icInd = 3; P.yrsInd = 3;

%     %These calculations won't work right if the invest cost interval isn't infinite
%     initBoreBen = P.idBoreInt*P.boreMax - P.idBoreSlope/2*P.boreMax.^2 - P.electricityBore.*P.initialLift.*P.boreMax;
%     initDugBen = P.idDugInt*P.dugMax - P.idDugSlope/2*P.dugMax.^2 - P.electricityDug.*P.initialLift.*P.dugMax;
% 
% %     P.investCostBase = (P.discount*(initBoreBen - initDugBen)-P.investCostNow*(1-P.discount))/(P.discount*P.icDecayRate); %set the common investment to facilitate the current share being optimal
% %     P.costIn = P.investCostNow - P.investCostBase; %the individual cost component for the last farm who (priv) optimally invested last period
% % 
% %     P.investCostMean = P.costIn - P.investCostSD*norminv(P.shrBore0);
% % 
% %     if P.costIn<0
% %         warning('Individual investment cost shocks are currently negative, might imply a parameter problem')
% %     end
% 
    P.inTruncProb = normcdf(P.maxInvestCost,P.investCostMean,P.investCostSD) - normcdf(P.minInvestCost,P.investCostMean,P.investCostSD);
    P.probBelow = normcdf(P.minInvestCost,P.investCostMean,P.investCostSD);
    P.lowCost = norminv(P.shrBore0*P.inTruncProb + P.probBelow,P.investCostMean,P.investCostSD);

%     if P.inTruncProb<1 || P.probBelow>0
%         error('The optimal in calcs assume no truncation')
%     end
    
    output.P = P;
end

%solve the rational expectations problem
 output.reOut = reSolve(P,modelOpts,P.h0*ones(modelOpts.minT,1));
%% solve the optimal management problem
model.func = 'optFunc';
model.discount = P.discount;
model.params = {P};
%model.horizon = 25;

smin(P.levelInd) = min(output.reOut.statePath(:,P.levelInd));
%smin(P.levelInd) = 10;
smax(P.levelInd) = P.landHeight;
smin(P.sbInd) = min(P.maxShr-.01,P.shrBore0); %allows the general code to work when shrBore is fixed at 1.
smax(P.sbInd) = P.maxShr;

smin(P.icInd) = 0;
smax(P.icInd) = P.investCostBase;

n = [modelOpts.capNodes modelOpts.heightNodes modelOpts.icNodes];
fspace = fundefn('lin',n,smin,smax);
minHeightDug = P.landHeight - P.maxDepthDug;

%Make sure there are nodes at maxDepthDug and the point where limits kick in dugWells
if P.landHeight-smin(P.levelInd)>P.maxDepthDug
    heightBreaks = fspace.parms{P.levelInd}{1};
    fspace.parms{P.levelInd}{1} = unique(sort([heightBreaks; P.maxDepthDug; P.landHeight-P.liftFullD]));
    fspace.parms{P.levelInd}{2} = 0;
    fspace.n(P.levelInd) = length(fspace.parms{P.levelInd}{1});
end
% 
% %Add a dense grid of points near begining and end of investment intervals
% capBreaks = fspace.parms{P.sbInd}{1};
% fspace.parms{P.sbInd}{1} = unique(sort([capBreaks; (P.shrBore0+.01:.01:P.shrBore0+.2)'; (.95:.01:.99)']));
% fspace.parms{P.sbInd}{2} = 0;
% fspace.n(P.sbInd) = length(fspace.parms{P.sbInd}{1});

snodes = funnode(fspace);
s = gridmake(snodes);

if exist(['solGuess' solveType '.mat'],'file')
    disp('loading old guess for opt control')
    load(['solGuess' solveType])
    vGuess = funeval(vGuessC,fspaceGuess,s);
    for ii=1:size(xGuessC,2)
        xGuess(:,ii) = funeval(xGuessC(:,ii),fspaceGuess,s);
    end
    if strfind(solveType,'niP')
        xGuess(:,P.investInd) = 0;
    end
else
    disp('generating new guess for opt control')
    %make guess at the optimal actions
    lifts = P.landHeight - s(:,P.levelInd);
    costBore = P.electricityBore*lifts;
    costDug = P.electricityDug*lifts;

    states = length(s);
    [lb,ub] = feval(model.func,'b',s,ones(states,3),[],P);

    xGuess(:,P.gwDugInd) = max((P.idDugInt - costDug)./P.idDugSlope,0);
    xGuess(:,P.gwBoreInd) = max((P.idBoreInt - costBore)./P.idBoreSlope,0);

    nb = netBen(xGuess(:,P.gwDugInd),xGuess(:,P.gwBoreInd),xGuess(:,P.investInd),s(:,P.levelInd),s(:,P.sbInd),s(:,P.icInd),P);
    deltaBen = nb.bore - nb.dug;
    newShr = max(s(:,P.sbInd),normcdf(deltaBen*P.discount/(1-P.discount),P.investCostMean,P.investCostSD));
    xGuess(:,P.investInd) = newShr-s(:,P.sbInd);

    xGuess = max(lb,min(xGuess,ub));
    vGuess = feval(model.func,'f',s,xGuess,[],P)/(1-P.discount);
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
if output.opt.converged == 1
    %save approximated fcns for value and action
    clear vGuessC xGuessC fspaceGuess
    vGuessC = c;
    fspaceGuess = fspace;
    dx = size(xGuess,2);
    xLong = reshape(x,numel(c),dx);
    for ii=1:dx
        xGuessC(:,ii) = funfitxy(fspace,funbasx(fspace,s),xLong(:,ii));
    end
    
    %save the solution to use as a starting point next time
    save(['solGuess' solveType],'vGuessC','xGuessC','fspaceGuess')
end

[shares,levels] = ndgrid(scoord{P.sbInd},scoord{P.levelInd});
s0(P.sbInd) = P.shrBore0;
s0(P.levelInd) = P.h0;
s0(P.icInd) = P.investCostBase;
[ssim,xsim] = dpsimul(model,s0,modelOpts.minT,scoord,x);
close

output.opt.vFunc = v;
output.opt.val = funeval(c,fspace,s0);
output.opt.statePath = squeeze(ssim)';
output.opt.controlPath = squeeze(xsim)';
optNb = netBen(output.opt.controlPath(:,P.gwDugInd),output.opt.controlPath(:,P.gwBoreInd),output.opt.controlPath(:,P.investInd),output.opt.statePath(:,P.levelInd),output.opt.statePath(:,P.sbInd),output.opt.statePath(:,P.icInd),P);
output.opt.valPath = [optNb.dug optNb.bore optNb.all];
output.opt.optVal = (P.discount.^(0:length(output.opt.valPath)-1))*output.opt.valPath(:,3);

output.reOut.pgain = (output.opt.val - output.reOut.reVal)/output.reOut.reVal;

%figure('Visible',modelOpts.figureVisible)
xTitles = {'Investment','(Traditional)','(Modern)'};
sTitles = {'Share in Modern Agriculture','Pumping Lift'};
sYlabel = {'%','meters'};
% 
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