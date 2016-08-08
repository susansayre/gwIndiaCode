function reOutput = reSolve(P,modelOpts,initialPath)

t = 1;

npvTol = 1e-3;
P.trendTol = 1e-1;
shrBore = P.shrBore0;

assumedPath = initialPath;
pathChange = P.trendTol*10;
iter = 0;
cycleCount = 1;
moveFrac = .1;
% optset('dpsolve','algorithm','funcit');
while pathChange>P.trendTol
    
    lift = P.landHeight - assumedPath;
    costDug = P.electricityDug*lift;
    costBore = P.electricityBore*lift;

    [lb,ub] = optFunc('b',[ones(size(assumedPath)) assumedPath],[],[],[],[],P);
    maxDugDemand = max(lb(:,P.gwDugInd),(P.idDugInt-costDug)./P.idDugSlope);
    gwDug = min(maxDugDemand,ub(:,P.gwDugInd));
    
    gwBore = min(max(lb(:,P.gwBoreInd),(P.idBoreInt-costBore)./P.idBoreSlope),ub(:,P.gwBoreInd));

    nbDug = P.idDugInt*gwDug - P.idDugSlope/2*gwDug.^2 - costDug.*gwDug;
    nbBore = P.idBoreInt*gwBore - P.idBoreSlope/2*gwBore.^2 - costBore.*gwBore;

    deltaBen = nbBore-nbDug;
    npvDeltaBen = (P.discount.^(0:1:length(deltaBen)-1)').*deltaBen;
    
    minPrefer2Yest = P.discount/(1-P.discount)*deltaBen(1:end-1);
    maxPrefer2Tom = P.discount/(1-P.discount)*deltaBen(2:end);
    
    for tt=1:length(assumedPath)-1;
        minPrefer21(tt,:) = 1/(1-P.discount^(tt-1))*sum(npvDeltaBen(2:tt));
        maxPrefer2N(tt,:) = 1/(P.discount^(tt-1))*sum(npvDeltaBen(tt+1:end));
    end
    
    minPrefer21(1) = 0; %guarantees that time 1 "prefers" itself to itself
    minCostInvestNow = max(minPrefer2Yest,minPrefer21);
    maxCostInvestNow = min(maxPrefer2Tom,maxPrefer2N);
  
    %keyboard
    levelPath = P.h0;
    shrPath = P.shrBore0;
    t=1;
    npvChange = 10*npvTol;
    investPath = [];
    valPath = [];
    while npvChange>npvTol
        lowCost = norminv(shrPath(t)*P.inTruncProb+P.probBelow,P.investCostMean,P.investCostSD);
        if t<length(gwBore)
           highCost = min(max(lowCost,maxCostInvestNow(t)),P.maxInvestCost);
           if t>1 && lowCost>minCostInvestNow(t)
             keyboard
           end
        else
           highCost = lowCost;
        end
        if t>=length(gwBore)
            gwBore(t+1) = gwBore(end);
            gwDug(t+1) = gwDug(end);
            nbDug(t+1) = nbDug(end);
            nbBore(t+1) = nbBore(end);
            assumedPath(t+1) = assumedPath(end);
        end
        gwTot = shrPath(t)*gwBore(t) + (1-shrPath(t))*gwDug(t);
        levelPath(t+1,:) = updateLevels(levelPath(t),gwTot,P);
        if highCost == lowCost;
            shrPath(t+1,:) = shrPath(t,:);
        else
            shrPath(t+1,:) = 1/P.inTruncProb*(normcdf(highCost,P.investCostMean,P.investCostSD)-P.probBelow);
        end
        investPath(t,:) = shrPath(t+1)-shrPath(t);
        if investPath(end)<0
            keyboard
        end
        investCost(t) = 1/P.inTruncProb*(P.investCostSD^2*(normpdf(lowCost,P.investCostMean,P.investCostSD)-normpdf(highCost,P.investCostMean,P.investCostSD))+P.investCostMean*investPath(t));
        valPath(t,:) = [nbDug(t) nbBore(t) shrPath(t)*nbBore(t) + (1-shrPath(t))*nbDug(t) - investCost(t)];
        npvChange = abs((P.discount^t)*valPath(t,3));
        t = t+1;
    end
    
%     figure()
%     plot(levelPath)
%     hold on;
%     plot(assumedPath,'--')
%     pause(1)
%     close
    pathChange = max(abs(assumedPath(1:length(levelPath))-levelPath));
    fprintf ('%4i %10.1e\n',iter,pathChange)
    
    cycleCount = cycleCount+1;
    if cycleCount>100
        moveFrac = moveFrac/2;
        cycleCount = 1;
    end
    assumedPath = (1-moveFrac)*assumedPath(1:length(levelPath)) +moveFrac*levelPath;
    iter = iter+1;
    %keyboard
    clear minPrefer21 maxPrefer2N
end

controlLength = length(investPath);
reOutput.statePath = [shrPath levelPath];
reOutput.controlPath = [investPath gwDug(1:controlLength) gwBore(1:controlLength)];
reOutput.valPath = valPath;
reOutput.reVal = (P.discount.^(0:length(valPath)-1))*valPath(:,3);
