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

assumedLevelPath = initialPath(1:end-1);

levelPathChange = trendTol*10;
iter = 0;

cycleCount = 1;
moveFrac = .1;
% optset('dpsolve','algorithm','funcit');
while levelPathChange>trendTol
    
    if iter>0
       clear maxPrefer2N maxPrefer2Tom minPrefer2Yest
    end
    
    lift = P.landHeight - assumedLevelPath;
    costDug = P.electricityDug*lift;
    costBore = P.electricityBore*lift;

    maxDugDemand = max(0,(P.idDugInt-costDug)./P.idDugSlope);
    dugUB = P.dugMax*(1-min(1,max(0,lift-P.liftFullD)/(P.maxDepthDug-P.liftFullD)));
    gwDug = min(maxDugDemand,dugUB);
    
    boreMax = P.boreMax - P.boreLimitDecline*max(0,lift-P.liftFullD); %ift lift is smaller than liftFullD, diff will be negative and we won't change limit
    gwBore = min(max(0,(P.idBoreInt-costBore)./P.idBoreSlope),boreMax);

    nbDug = P.idDugInt*gwDug - P.idDugSlope/2*gwDug.^2 - costDug.*gwDug;
    nbBore = P.idBoreInt*gwBore - P.idBoreSlope/2*gwBore.^2 - costBore.*gwBore;
    
    exponents = (1:length(gwBore))';
    dFs = P.discount.^exponents;
    deltaBen = nbBore-nbDug;
 
    deltaBenNext = [deltaBen(2:end); deltaBen(end)/(1-P.discount)];
    npvDeltaBenNext = (P.discount.^(exponents+1)).*deltaBenNext;
    maxPrefer2Tom = (P.discount*deltaBenNext - (1-P.discount*(1-P.icDecayRate))*P.investCostBase*((1-P.icDecayRate).^exponents))./(1-P.discount);
    
    minPrefer2Yest = (P.discount*deltaBen - (1-P.discount*(1-P.icDecayRate))*P.investCostBase*((1-P.icDecayRate).^(exponents-1)))./(1-P.discount);

    maxPrefer2N = Inf + zeros(length(deltaBen),1);
    for tt=1:length(deltaBen);
        maxPrefer2N(tt,:) = 1/(P.discount^tt)*sum(npvDeltaBenNext(tt:end))-((1-P.icDecayRate)^tt)*P.investCostBase;
    end

    maxCostInvestNow = min(maxPrefer2Tom,maxPrefer2N); %if my cost is above this, I should either invest tomorrow or never
    minCostInvestNow = minPrefer2Yest;
    doubleIntervals = 1;
    while doubleIntervals>0
        %identify intervals that are possibly optimal for investment for some
        %farms
        possibleInvestIntervals = find(minCostInvestNow<maxCostInvestNow);
        
        %check whether the high and mid cost values for each interval are also
        %in another interval (note that the low cost value for each interval is
        %the high cost value for the previous interval by definition)

        maxVals = maxCostInvestNow(possibleInvestIntervals);
        minVals = minCostInvestNow(possibleInvestIntervals);
        maxValsMat = repmat(maxVals,1,length(maxVals));
        minValsMat = repmat(minVals,1,length(maxVals));

        inIntervals = (maxValsMat'<=maxValsMat).*(maxValsMat'>minValsMat);
        problemCosts = find(sum(inIntervals)>1);
        doubleIntervals = numel(problemCosts);
        
        for ci=1:numel(problemCosts)
            %keyboard
            chkIntervals = find(inIntervals(:,problemCosts(ci)));
            time1 = possibleInvestIntervals(problemCosts(ci));
            maxPreferNow  = maxCostInvestNow(time1);     
            minPreferNow = minCostInvestNow(time1);
            for ti=1:numel(chkIntervals)
                time2 = possibleInvestIntervals(chkIntervals(ti));
                if time1>time2
                    minPreferNow = max(minPreferNow,(sum(dFs(time2+1:time1).*deltaBen(time2+1:time1)) - P.investCostBase*((P.discount*(1-P.icDecayRate))^time2 - (P.discount*(1-P.icDecayRate))^time1))/(P.discount^time2 - P.discount^time1));
                    maxCostInvestNow(time2) = min(maxCostInvestNow(time2),minPreferNow);
                elseif time1==time2
                    continue
                else
                    maxPreferNow = min(maxPreferNow,(sum(dFs(time1+1:time2).*deltaBen(time1+1:time2)) - P.investCostBase*((P.discount*(1-P.icDecayRate))^time1 - (P.discount*(1-P.icDecayRate))^time2))/(P.discount^time1 - P.discount^time2));
                    minCostInvestNow(time2) = maxPreferNow;
                end
            end
            maxCostInvestNow(time1)= maxPreferNow;
            minCostInvestNow(time1) = minPreferNow; 
        end
    end
            
    %simulate the consequences of the implied behavior above
    %keyboard
    levelPath = P.h0;
    shrPath = P.shrBore0;
    invCostC = P.investCostBase*(1-P.icDecayRate); %for consistency with later code
    t=1;
    npvChange = 10*npvTol;
    investPath = [];
    valPath = [];
    while npvChange>npvTol || t<modelOpts.minT%continues simulating forward until the NPV contribution is negligible
        lowCost = norminv(shrPath(t)*P.inTruncProb+P.probBelow,P.investCostMean,P.investCostSD);     
        if t<=length(gwBore)
           %the max operator will guarantee that we don't want to de-invest if I've already invested but wouldn't want to now
           %the min operator deals with the truncation in the investCost distribution
           if maxCostInvestNow(t)<=minCostInvestNow(t)
               highCost = lowCost;
           else
               highCost = min(max(lowCost,maxCostInvestNow(t)),P.maxInvestCost);
           end
%            if t>1 && lowCost>minCostInvestNow(t)
%              keyboard
%            end
            %recompute use based on current depth
            thisLift = P.landHeight - levelPath(t);
            thisCostDug = P.electricityDug*thisLift;
            thisCostBore = P.electricityBore*thisLift;

            thisMaxDugDemand = max(0,(P.idDugInt-thisCostDug)./P.idDugSlope);
            thisDugUB = P.dugMax*(1-min(1,max(0,thisLift-P.liftFullD)/(P.maxDepthDug-P.liftFullD)));
            thisGwDug = min(thisMaxDugDemand,thisDugUB);
    
            thisBoreMax = P.boreMax - P.boreLimitDecline*max(0,thisLift-P.liftFullD); %ift lift is smaller than liftFullD, diff will be negative and we won't change limit
            thisGwBore = min(max(0,(P.idBoreInt-thisCostBore)./P.idBoreSlope),thisBoreMax);

            thisNbDug = P.idDugInt*thisGwDug - P.idDugSlope/2*thisGwDug.^2 - thisCostDug.*thisGwDug;
            thisNbBore = P.idBoreInt*thisGwBore - P.idBoreSlope/2*thisGwBore.^2 - thisCostBore.*thisGwBore;

       else
            %I've hit the end of my data, keep levels (and therefore net benefits) constant at their last levels and
            %prevent investment going forward, the next iteration will proceed this far if we need it and see how
            %people respond to the actual levels induced by this behavior
            gwBore(t) = gwBore(end);
            gwDug(t) = gwDug(end);
            nbDug(t) = nbDug(end);
            nbBore(t) = nbBore(end);
            assumedLevelPath(t) = assumedLevelPath(end);
            highCost = lowCost;
        end
        gwTot = shrPath(t)*thisGwBore + (1-shrPath(t))*thisGwDug;
        levelPath(t+1,:) = updateLevels(levelPath(t),gwTot,P);
        invCostC(t+1,:) = invCostC(t)*(1-P.icDecayRate);
        if highCost == lowCost; %doing it this way prevents weird rounding problems
            shrPath(t+1,:) = shrPath(t,:);
        else
            shrPath(t+1,:) = 1/P.inTruncProb*(normcdf(highCost,P.investCostMean,P.investCostSD)-P.probBelow);
        end
        investPath(t,:) = shrPath(t+1)-shrPath(t);
        if investPath(t,:)>P.investLimit
            investPath(t,:) = P.investLimit;
            shrPath(t+1,:) = shrPath(t,:) + investPath(t,:);
            highCost = norminv(shrPath(t+1,:)*P.inTruncProb+P.probBelow,P.investCostMean,P.investCostSD);
        end
        if investPath(end)<0
            keyboard
        end
        investCost(t) = 1/P.inTruncProb*(P.investCostSD^2*(normpdf(lowCost,P.investCostMean,P.investCostSD)-normpdf(highCost,P.investCostMean,P.investCostSD))+P.investCostMean*investPath(t)) + investPath(t)*invCostC(t);
        valPath(t,:) = [thisNbDug thisNbBore shrPath(t)*thisNbBore + (1-shrPath(t))*thisNbDug - investCost(t)];
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
%     figure()
%     plot(levelPath)
%     hold on;
%     plot(assumedLevelPath,'--')
%     pause(.001)
%     close
    
    assumedLevelPath = (1-moveFrac)*assumedLevelPath(1:t-1) + moveFrac*levelPath(1:t-1); 
    iter = iter+1;
    %keyboard

end

%check optimal path is globally, not locally optimal
investIntervals = find(investPath>0);
lowCost = norminv(shrPath(investIntervals),P.investCostMean,P.investCostSD);
highCost = norminv(shrPath(investIntervals+1),P.investCostMean,P.investCostSD);
midCost = (lowCost + highCost)/2;
dFs = P.discount.^(1:length(nbDug))';

for ti=1:t-1
    npvInvestCommon(ti,:) = sum(dFs(1:ti).*nbDug(1:ti)) - dFs(ti)*invCostC(ti) + sum(dFs(ti+1:end).*nbBore(ti+1:end)) + dFs(end)*P.discount/(1-P.discount)*nbBore(end);
end
npvInvestCommon(t,:) = sum(dFs.*nbDug) + dFs(end)*nbDug(end)*P.discount/(1-P.discount);
dFHat = [dFs(1:t-1); 0];
for jj=1:numel(investIntervals);
    npvInvestLow(:,jj) = npvInvestCommon - dFHat*lowCost(jj);
    npvInvestHigh(:,jj) = npvInvestCommon - dFHat*highCost(jj);
    npvInvestMid(:,jj) = npvInvestCommon - dFHat*midCost(jj);
end

if any(investIntervals)
    [maxValMid,maxYrMid] = max(npvInvestMid);
    [maxValHigh,maxYrHigh] = max(npvInvestHigh);
    [maxValLow,maxYrLow] = max(npvInvestLow);

    if any(maxYrMid-investIntervals')
        keyboard
    end
end
controlLength = length(investPath);
reOutput.statePath = [shrPath levelPath];
reOutput.controlPath = [investPath gwDug(1:controlLength) gwBore(1:controlLength)];
reOutput.valPath = valPath;
reOutput.reVal = (P.discount.^(0:length(valPath)-1))*valPath(:,3);

