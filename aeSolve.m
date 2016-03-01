function aeOutput = aeSolve2(P,modelOpts)
%solve the adaptive expectations problem as an optimal switching problem

statePath(1,P.levelInd) = P.h0; 
statePath(1,P.sbInd) = P.shrBore
trendPath = P.levelTrend

t=1;
change = 1;
P.itol = 1e-5;

while t<P.minT||change>P.iTol
    
%compute the optimal action this year and next year for each type given the
%trajectory.
    levels = [statePath(t,P.levelInd); statePath(t,P.levelInd)-trendPath(t)];
    shrBore = statePath(t,P.sbInd);
    lift = P.landHeight - levels;
    costDug = P.costDug_a*exp(P.costDug_b*lift);
    costBore = P.electricity*lift;

    gwDug = (P.dDugInt-costDug)/P.dDugSlope;
    gwBore = (P.dBoreInt-costBore)/P.dDugBore;
    controlPath(t,P.gwDugInd) = gwDug;
    controlPath(t,P.gwBoreInd) = gwBore;
    
    nbDug = P.dDugInt*gwDug - P.dDugSlope/2*gwDug.^2 - costDug.*gwDug;
    nbBore = P.dBoreInt*gwBore - P.dBoreSlope/2*gwBore.^2 - costBore.*gwBore;
    
    deltaBen = nbBore(2) - nbDug(2);
    newShr = normcdf(P.discount*deltaBen,P.investCostMean,P.investCostSD);
    invest = newShr - shrBore;
    controlPath(t,P.investInd) = invest;
    
    lowCost = min(max(norminv(shrBore,P.investCostMean,P.investCostSD),-1/eps),1/eps);
    highCost = min(max(norminv(newShr,P.investCostMean,P.investCostSD),-1/eps),1/eps);
    deltaCost = highCost-lowCost;
    
%     approxShare = [0:.1:1]
%     costArray = repmat(lowCost,1,size(approxShare,2)) + repmat(deltaCost,1,size(approxShare,2)).*approxShare;
%     probArray = normpdf(costArray,P.investCostMean,P.investCostSD);
%     intervalIntegral = sum(costArray.*probArray,2);
    
    intervalIntegral = integral(@(x) x.*normpdf(x,P.investCostMean,P.investCostSD),lowCost,highCost);
    
    investCost = intervalIntegral.*invest; 
    
    value(t,:) = [nbDug(1) nbBore(1) shrBore*nbBore(1)+(1-shrBore)*nbDug(1)-investCost];
    npvValue(t,:) = value(t,:)*P.discount^(t-1);

    statePath(t+1,P.levelInd) = updateLevels(level(1),shrBore*gwBore+(1-shrBore)*gwDug,P);
    statePath(t+1,P.sbInd) = newShr;
    
    change = max(invest,abs(statePath(t+1,P.levelInd)-levels(2)));
    t = t+1;
end

aeOutput.statePath = statePath;
aeOutput.controlPath = controlPath;
aeOutput.valPath = value;
aeOutput.value = sum(npvValue);
   
	
