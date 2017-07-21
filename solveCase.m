function [output] = solveCase(P,modelOpts,parameterSetID)

if strfind(parameterSetID{2},'only2')
    solveType = 'only2';
elseif strfind(parameterSetID{2},'only1')
    solveType = 'only1';
else
    solveType = '';
end

myFields = fieldnames(P);
for ii=1:numel(myFields)
    thisParam = myFields{ii};
    if strcmp(thisParam,'csString')||strcmp(thisParam,'viableTechs')
        continue
    end
    if eval(['numel(P.' thisParam ')>numel(P.viableTechs)'])
        eval(['P.' thisParam '= P.' thisParam '(P.viableTechs);'])
    end
end
P.numTech = numel(P.viableTechs);
P.idInts = P.chokePrice*ones(1,P.numTech);
slope = (P.chokePrice - P.mbAtMax)/P.maxWater;
P.idSlopes = slope*ones(1,P.numTech);
P.numActionOpt = P.numTech^2;
%create indices
P.ind.level = 1;
P.shareInds = 2:P.numTech;

P.ind.icM = 1;
P.ind.tech = 2;
P.ind.yrInd = 3;

P.ind.water = 1:P.numTech;
P.ind.invest = P.numTech + (1:(P.numTech*(P.numTech-1)));

P.interactive = 0;

%compute derived parameter values
P.eCosts = (P.energyNeed./P.pumpEfficiency).*P.kWhCost;

P.fixedEcosts = (P.energyNeed./P.pumpEfficiency).*P.maxLiftMaxCap*P.maxWater*P.avgPerKwH;

P.h0 = P.landHeight - P.initialLift;
P.inflow = P.maxWater*(1-P.returnFlow)/P.maxExploitRatio;
P.AS = (P.currIrrigShare*P.maxWater*(1-P.returnFlow) - P.inflow)/P.annualDrop;
P.s0(:,P.ind.level) = P.h0;
P.s0(:,P.shareInds) = share2RawShr(P.initShares);

P.maxICM = P.maxIC/P.baseIC - 1;
P.maxProfit = (0.5*P.chokePrice+0.5*P.mbAtMax)*P.maxWater;
P.baseICyears = P.baseIC/P.maxProfit;
P.maxICyears = P.maxIC/P.maxProfit;
for ii=1:P.numTech
    P.techCostMat(ii,:) = eval(['P.baseIC*P.ic' num2str(ii) 'cost']);
end
%keyboard
P.transitionMat = repmat(1:P.numTech,P.numTech,1);
[P.orderedTechCosts,P.techCostOrder] = sort(P.techCostMat,2,'descend');
%solve the rational expectations problem

if modelOpts.refineExisting
    load(['detailedOutput/' parameterSetID{1} '/' parameterSetID{2} ''],'output')
    initialPath = output.reCp.levelPath;
    prevPathReRealE = output.reRealE.levelPath;
    clear output
else
    initialPath = P.h0-(0:modelOpts.minT)';
end
output.reCp = reSolve(P,modelOpts,initialPath);

PrealE = P;
PrealE.eCostShr = 1;
PrealE.fixedEcosts = 0*P.fixedEcosts;

if modelOpts.refineExisting
    initialPath = prevPathReRealE;
end
output.reRealE = reSolve(PrealE,modelOpts,initialPath);

% PbigT = P;
% PbigT.fixedEcosts = 4*P.fixedEcosts;
% output.reBigT = reSolve(PbigT,modelOpts,initialPath);
% 
% PsubE = PrealE;
% PsubE.eCostShr = sum(output.reCp.nbDetail(:,7))/sum(output.reCp.nbDetail(:,2));
% output.reSubE = reSolve(PsubE,modelOpts,initialPath);

save reCp
%% solve the optimal management problem
model.func = 'optFunc';
model.discount = P.discount;
model.params = {P};
%model.horizon = 25;

smin(P.ind.level) = min(output.reRealE.levelPath);
smax(P.ind.level) = P.landHeight-P.minLift;
smin(P.shareInds) = 0;
smax(P.shareInds) = 1;

n(P.ind.level) = modelOpts.heightNodes;
n(P.shareInds) = modelOpts.capNodes;

fspace = fundefn('lin',n,smin,smax);

%Make sure there are nodes at maxDepthDug and the point where limits kick in dugWells
for jj=1:P.numTech

    if P.landHeight-smin(P.ind.level)>P.maxDepths(jj)
        heightBreaks = fspace.parms{P.ind.level}{1};
        fspace.parms{P.ind.level}{1} = unique(sort([heightBreaks; P.landHeight - P.maxDepths(jj)]));
        fspace.parms{P.ind.level}{2} = 0;
        fspace.n(P.ind.level) = length(fspace.parms{P.ind.level}{1});
    end
end

snodes = funnode(fspace);
s = gridmake(snodes);
[ns,ds] = size(s);

if modelOpts.refineExisting
    disp('loading previous solution')
    load(['detailedOutput/' parameterSetID{1} '/' parameterSetID{2} ''],'vGuessC','xGuessC','fspaceGuess')
    vGuess = funeval(vGuessC,fspaceGuess,s);
    for ii=1:size(xGuessC,2)
        xGuess(:,ii) = funeval(xGuessC(:,ii),fspaceGuess,s);
    end
elseif exist(['solGuess' solveType '.mat'],'file')
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
    lifts = P.landHeight - s(:,P.ind.level); %ns x 1
    costs = lifts*P.eCosts;

    states = length(s);
    [lb,ub] = feval(model.func,'b',s,ones(states,3),[],P);

    Intercepts = repmat(P.idInts,ns,1);
    Slopes = repmat(P.idSlopes,ns,1);

    waterGuess = max((Intercepts - costs)./Slopes,0);

    %set investGuess so that 90% of the farms on any viable type remain in that type
    % 90% of the farms that convert only go up one tech, and so on
    % no farms go backward and all farms upgrade at least one tech when a technology is non-viable
    investGuess = zeros(ns,P.numTech*(P.numTech-1));
    investGuessHat = investGuess;
    investGuess(find(s(:,P.shareInds(1))),1:P.numTech-1) = .9;
    investGuessHat(find(s(:,P.shareInds(1))),1) = 1;
    for ii=2:P.numTech-1
        investGuess(find(s(:,P.shareInds(ii))),(ii-1).*(P.numTech-1)+ii:ii*(P.numTech-1)) = .9;
        investGuessHat(find(s(:,P.shareInds(ii))),(ii-1).*(P.numTech-1)+ii) = 1;
    end
    xGuess(:,P.ind.water) = waterGuess;
    xGuess(:,P.ind.invest) = investGuess;
    xGuessHat = xGuess;
    xGuessHat(:,P.ind.invest) = investGuessHat;
    vGuess = (1-P.discount)/P.discount*optFunc('f',s,xGuessHat,[],P);
end

optset('dpsolve','algorithm',modelOpts.algorithm);
optset('dpsolve','maxit',modelOpts.maxit);
optset('dpsolve','showiters',1);
optset('dpsolve','tol',modelOpts.vtol);
optset('dpsolve','nres',modelOpts.nres);
optset('dpsolve_vmax','maxbacksteps',0);
optset('dpsolve_vmax','maxit',50);
optset('dpsolve_vmax','lcpmethod','ss');

%figure()
save beforeOpts
[c,scoord,v,x,resid,exf] = dpsolveS2(model,fspace,s,vGuess,xGuess);
%[c,scoord,v,x] = dpsolve(model,fspace,s,vGuess,xGuess);

output.opt.converged = exf;
if output.opt.converged == 1
    %save approximated fcns for value and action
    clear vGuessC xGuessC fspaceGuess
    vGuessC = c;
    fspaceGuess = fspace;
    sVals = gridmake(scoord);
    dx = size(xGuess,2);
    xLong = reshape(x,size(sVals,1),dx);
    for ii=1:dx
        xGuessC(:,ii) = funfitxy(fspace,sVals,xLong(:,ii));
    end
    
   %save the solution to use as a starting point next time
   save(['solGuess' solveType],'vGuessC','xGuessC','fspaceGuess')
end

save afterOpt

s0 = P.s0;
[ssim,xsim] = dpsimul(model,s0,modelOpts.minT-1,scoord,x);

output.opt.vFunc = v;
output.opt.val = funeval(c,fspace,s0);
output.opt.statePath = squeeze(ssim)';
output.opt.sharePath = rawShr2Share(output.opt.statePath(:,P.shareInds));
output.opt.levelPath = output.opt.statePath(:,P.ind.level);
output.opt.controlPath = squeeze(xsim)';
[nbS,nbDetailS,fDetailS] = optValue(output.opt.statePath,output.opt.controlPath,[],P,'society');
[nbI,nbDetailI,fDetailsI] = optValue(output.opt.statePath,output.opt.controlPath,[],P,'indiv');
output.opt.nbPath = [nbS nbI];
output.opt.nbDetailPath = [nbDetailS nbDetailI];
output.opt.npv = P.discount.^(0:length(output.opt.nbPath)-1)*output.opt.nbPath;

if P.checkAtInvest && P.maxInvest
    %identify points when investment threshold is reached
    realEYear = min(find(1-output.reRealE.sharePath(:,1)>=P.checkAtInvest));
    if realEYear
        s0RealETest = output.reRealE.statePath(realEYear,:);
        [ssimRealETest,xsimRealETest] = dpsimul(model,s0RealETest,modelOpts.minT-realEYear-1,scoord,x);
        output.opt.realEStartTest.statePath = squeeze(ssimRealETest)';
        output.opt.realEStartTest.sharePath = rawShr2Share(output.opt.realEStartTest.statePath(:,P.shareInds));
        output.opt.realEStartTest.levelPath = output.opt.realEStartTest.statePath(:,P.ind.level);
        output.opt.realEStartTest.controlPath = squeeze(xsimRealETest)';
        [nbS,nbDetailS] = optValue(output.opt.realEStartTest.statePath,output.opt.realEStartTest.controlPath,[],P,'society');
        [nbI,nbDetailI] = optValue(output.opt.realEStartTest.statePath,output.opt.realEStartTest.controlPath,[],P,'indiv');
        output.opt.realEStartTest.nbPath = [nbS nbI];
        output.opt.realEStartTest.nbDetailPath = [nbDetailS nbDetailI];
        output.opt.realEStartTest.npv = P.discount.^(0:length(output.opt.realEStartTest.nbPath)-1)*output.opt.realEStartTest.nbPath;

        output.reRealE.realEStartTest.npv = P.discount.^(0:length(output.opt.realEStartTest.nbPath))*output.reRealE.nbPath(realEYear:end,:);
        output.reRealE.realEStartTest.pgain  = (output.opt.realEStartTest.npv - output.reRealE.realEStartTest.npv)./output.reRealE.realEStartTest.npv;
    end

    reCpYear = min(find(1-output.reCp.sharePath(:,1)>=P.checkAtInvest));
    if reCpYear
        s0CpTest = output.reCp.statePath(reCpYear,:);
        [ssimCpTest,xsimCpTest] = dpsimul(model,s0CpTest,modelOpts.minT-reCpYear-1,scoord,x);
        output.opt.reCpStartTest.statePath = squeeze(ssimCpTest)';
        output.opt.reCpStartTest.sharePath = rawShr2Share(output.opt.reCpStartTest.statePath(:,P.shareInds));
        output.opt.reCpStartTest.levelPath = output.opt.reCpStartTest.statePath(:,P.ind.level);
        output.opt.reCpStartTest.controlPath = squeeze(xsimCpTest)';
        [nbS,nbDetailS] = optValue(output.opt.reCpStartTest.statePath,output.opt.reCpStartTest.controlPath,[],P,'society');
        [nbI,nbDetailI] = optValue(output.opt.reCpStartTest.statePath,output.opt.reCpStartTest.controlPath,[],P,'indiv');
        output.opt.reCpStartTest.nbPath = [nbS nbI];
        output.opt.reCpStartTest.nbDetailPath = [nbDetailS nbDetailI];
        output.opt.reCpStartTest.npv = P.discount.^(0:length(output.opt.reCpStartTest.nbPath)-1)*output.opt.reCpStartTest.nbPath;

        output.reCp.reCpStartTest.npv = P.discount.^(0:length(output.opt.reCpStartTest.nbPath))*output.reCp.nbPath(reCpYear:end,:);
        output.reCp.reCpStartTest.pgain  = (output.opt.reCpStartTest.npv - output.reCp.reCpStartTest.npv)./output.reCp.reCpStartTest.npv;
    end

end
    
reRuns = fieldnames(output);
for ii=1:numel(reRuns)
    if strcmp(reRuns{ii},'opt')
       continue
    else
        eval(['output.' reRuns{ii} '.pgain = (output.opt.npv - output.' reRuns{ii} '.npv)./output.' reRuns{ii} '.npv;'])
    end    
end

%% store details for later reference if needed
if ~exist('detailedOutput','dir')
    mkdir('detailedOutput')
end
if ~exist(['detailedOutput/' parameterSetID{1}],'dir')
    mkdir(['detailedOutput/' parameterSetID{1}])
end

yrs = modelOpts.minT-1;
plotResults

%plot re policy function
figure()
output.reCp.choicePlot.draw();
saveas(gcf,fullfile('detailedOutput',parameterSetID{1},[parameterSetID{2} '_reCpChoices']),'epsc')

figure()
output.reRealE.choicePlot.draw();
saveas(gcf,fullfile('detailedOutput',parameterSetID{1},[parameterSetID{2} '_reRealEChoices']),'epsc')
close all

save(['detailedOutput/' parameterSetID{1} '/' parameterSetID{2} ''])