function reOutput = reSolve(P,modelOpts,initialPath)

% model.func = 'optStoppingFunc';
% model.discount = P.discount;
% model.actions = [0;1];
% model.discretestates = 3;
% model.e = 0;
% model.w = 1;

npvTol = modelOpts.vtol;
trendTol = modelOpts.ttol;
shrTol = modelOpts.stol;

assumedLevelPath = initialPath(1:end-1,P.levelInd);
assumedShrPath = initialPath(1:end-1,P.sbInd);
assumedInvestPath = assumedShrPath(2:end) - assumedShrPath(1:end-1);

levelPathChange = trendTol*10;

iter = 0;

cycleCount = 1;
moveFrac = .1;
% optset('dpsolve','algorithm','funcit');
while levelPathChange>trendTol
    
    lift = P.landHeight - assumedLevelPath;
    costDug = P.electricityDug*lift;
    costBore = P.electricityBore*lift;

    maxDugDemand = max(0,(P.idDugInt-costDug)./P.idDugSlope);
    dugUB = P.dugMax*(1-min(1,max(0,lift-P.liftFullD)/(P.maxDepthDug-P.liftFullD)));
    gwDug = min(maxDugDemand,dugUB);
    gwBore = min(max(0,(P.idBoreInt-costBore)./P.idBoreSlope),P.boreMax);

    nbDug = P.idDugInt*gwDug - P.idDugSlope/2*gwDug.^2 - costDug.*gwDug;
    nbBore = P.idBoreInt*gwBore - P.idBoreSlope/2*gwBore.^2 - costBore.*gwBore;

    deltaBen = nbBore-nbDug;
    npvDeltaBen = (P.discount.^(1:length(deltaBen))').*deltaBen; %same length as assumed paths

    maxPrefer2Tom = (P.discount*deltaBen - (1-P.discount)*P.investCostBase)./((1-P.discount*(1-P.icDecayRate))*(1-P.icDecayRate).^((1:length(deltaBen))'));

    maxPrefer2N(length(deltaBen),:) = 1/((P.discount*(1-P.icDecayRate))^length(deltaBen))*(npvDeltaBen(end)*P.discount/(1-P.discount)-P.investCostBase);
    for tt=1:length(deltaBen)-1;
        maxPrefer2N(tt,:) = 1/((P.discount*(1-P.icDecayRate))^tt)*(sum(npvDeltaBen(tt+1:end)) + npvDeltaBen(end)*P.discount/(1-P.discount)-P.investCostBase);
    end

    maxCostInvestNow = min(maxPrefer2Tom,maxPrefer2N); %if my cost is above this, I should either invest tomorrow or never

    %simulate the consequences of the implied behavior above
    %keyboard
    levelPath = P.h0;
    shrPath = P.shrBore0;
    t=1;
    npvChange = 10*npvTol;
    investPath = [];
    valPath = [];
    while npvChange>npvTol || t<modelOpts.minT%continues simulating forward until the NPV contribution is negligible
        lowCost = norminv(shrPath(t)*P.inTruncProb+P.probBelow,P.investCostMean,P.investCostSD);     
        if t<=length(gwBore)
           %the max operator will guarantee that we don't want to de-invest if I've already invested but wouldn't want to now
           %the min operator deals with the truncation in the investCost distribution
           highCost = min(max(lowCost,maxCostInvestNow(t)),P.maxInvestCost); 
%            if t>1 && lowCost>minCostInvestNow(t)
%              keyboard
%            end
       else
            %I've hit the end of my data, keep levels (and therefore net benefits) constant at their last levels and
            %prevent investment going forward, the next iteration will proceed this far if we need it and see how
            %people respond to the actual levels induced by this behavior
            gwBore(t) = gwBore(end);
            gwDug(t) = gwDug(end);
            nbDug(t) = nbDug(end);
            nbBore(t) = nbBore(end);
            assumedLevelPath(t) = assumedLevelPath(end);
            assumedInvestPath(t) = 0;
            highCost = lowCost;
        end
        gwTot = shrPath(t)*gwBore(t) + (1-shrPath(t))*gwDug(t);
        levelPath(t+1,:) = updateLevels(levelPath(t),gwTot,P);
        if highCost == lowCost; %doing it this way prevents weird rounding problems
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
    
    levelPathChange = max(abs(assumedLevelPath(1:t-1)-levelPath(1:t-1)));
    cycleCount = cycleCount+1;
    if cycleCount>100
        moveFrac = moveFrac/2;
        cycleCount = 1;
    end

    fprintf ('%4i %10.1e\n',iter,levelPathChange)
    figure()
    subplot(1,2,1)
    plot(levelPath)
    hold on;
    plot(assumedLevelPath,'--')
    subplot(1,2,2)
    plot(shrPath)
    hold on;
    plot(assumedShrPath,'--')
    pause(1)
    close
    
    assumedLevelPath = (1-moveFrac)*assumedLevelPath(1:t-1) +moveFrac*levelPath(1:t-1);
    
    iter = iter+1;
    %keyboard
    clear maxPrefer2N maxPrefer2Tom
end

controlLength = length(investPath);
reOutput.statePath = [shrPath levelPath];
reOutput.controlPath = [investPath gwDug(1:controlLength) gwBore(1:controlLength)];
reOutput.valPath = valPath;
reOutput.reVal = (P.discount.^(0:length(valPath)-1))*valPath(:,3);
