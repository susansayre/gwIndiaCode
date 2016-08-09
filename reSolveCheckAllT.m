function reOutput = reSolve(P,modelOpts,initialPath)

% t = 1;
% model.func = 'optStoppingFunc';
% model.discount = P.discount;
% model.actions = [0;1];
% model.discretestates = 3;
% model.e = 0;
% model.w = 1;

npvTol = 1e-3;
P.trendTol = 1e-1;

assumedLevelPath = initialPath(:,P.levelInd);
assumedShrPath = initialPath(:,P.sbInd);
assumedInvestPath = assumedShrPath(2:end) - assumedShrPath(1:end-1);
    
investCostNodes = P.shrBore0:(1-P.shrBore0)/(modelOpts.capNodes-1):1;
investCostDist = norminv(P.inTruncProb*investCostNodes+P.probBelow,P.investCostMean,P.investCostSD); %possible investment costs of farms

pathChange = P.trendTol*10;
iter = 0;
cycleCount = 1;
moveFrac = .05;
% optset('dpsolve','algorithm','funcit');
while pathChange>P.trendTol
    
    lift = P.landHeight - assumedLevelPath;
    costDug = P.electricityDug*lift;
    costBore = P.electricityBore*lift;

    maxDugDemand = max(0,(P.idDugInt-costDug)./P.idDugSlope);
    dugUB = P.dugMax*(1-min(1,max(0,lift-P.depthFullD)));
    gwDug = min(maxDugDemand,dugUB);
    gwBore = min(max(0,(P.idBoreInt-costBore)./P.idBoreSlope),P.boreMax);

    nbDug = P.idDugInt*gwDug - P.idDugSlope/2*gwDug.^2 - costDug.*gwDug;
    nbBore = P.idBoreInt*gwBore - P.idBoreSlope/2*gwBore.^2 - costBore.*gwBore;
    
    discountFactors = (P.discount.^(0:1:length(nbDug)-1)'); 
    deltaBen = nbBore-nbDug;
    npvDeltaBen = (P.discount.^(0:1:length(deltaBen)-1)').*deltaBen;
 
    currentlyTrad = tril(ones(numel(assumedLevelPath))); %lower triangular matrix
    benConvertT = currentlyTrad*(nbDug.*discountFactors) + (1-currentlyTrad)*(nbBore.*discountFactors);

    investCostAdder = P.investCostPenalty*assumedInvestPath.^2;
    investCostAdder(end+1) = 0; %final row is the never convert option

    benLessInvestAdder = benConvertT - discountFactors.*investCostAdder;
    
    [netBen,optConvertTime] = max(repmat(benLessInvestAdder,1,numel(investCostDist)) - discountFactors*investCostDist);
  
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
           parcelsConvertingNow = find(optConvertTime==t);
           if any(parcelsConvertingNow)
               highCost = min(max(lowCost,max(investCostDist(parcelsConvertingNow))),P.maxInvestCost);
           else
               highCost = lowCost;
           end
        else
           highCost = lowCost;
        end
        if t>=length(gwBore)
            gwBore(t+1) = gwBore(end);
            gwDug(t+1) = gwDug(end);
            nbDug(t+1) = nbDug(end);
            nbBore(t+1) = nbBore(end);
            assumedLevelPath(t+1) = assumedLevelPath(end);
            assumedInvestPath(t+1) = assumedInvestPath(end);
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
    
    figure()
    subplot(1,2,1)
    plot(levelPath)
    hold on;
    plot(assumedLevelPath,'--')
    subplot(1,2,2)
    plot(investPath)
    hold on;
    plot(assumedInvestPath,'--')
    pause(1)
    close
    pathChange = max(abs(assumedLevelPath(1:length(levelPath))-levelPath)) + max(abs(assumedInvestPath(1:length(investPath))-investPath));
    fprintf ('%4i %10.1e\n',iter,pathChange)
    
    cycleCount = cycleCount+1;
    if cycleCount>100
        moveFrac = moveFrac/2;
        cycleCount = 1;
    end
    assumedLevelPath = (1-moveFrac)*assumedLevelPath(1:length(levelPath)) +moveFrac*levelPath;
    assumedInvestPath = (1-moveFrac)*assumedInvestPath(1:length(investPath)) + moveFrac*investPath;
    iter = iter+1;
    %keyboard
    clear minPrefer21 maxPrefer2N
end

controlLength = length(investPath);
reOutput.statePath = [shrPath levelPath];
reOutput.controlPath = [investPath gwDug(1:controlLength) gwBore(1:controlLength)];
reOutput.valPath = valPath;
reOutput.reVal = (P.discount.^(0:length(valPath)-1))*valPath(:,3);
