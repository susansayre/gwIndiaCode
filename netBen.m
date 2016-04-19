function [nb,dnb,ddnb] = netBen(gwDug,gwBore,invest,gwLevel,shrBore,P)
    %this function returns the aggregate net benefit as a function of gw
    %extraction (which now needs to be a two dimensional vector),
    %investment in regime switching, gwLevels, and the current share of
    %plots that are irrigated by bore wells;
    
    ns = length(shrBore);
    
    gwLevel = max(P.bottom,min(P.landHeight,gwLevel));
	lift = P.landHeight - gwLevel;
    
    if any(lift>P.landHeight-P.bottom)||any(lift<0)
        keyboard
    end
    
    %compute cost per unit for dug wells
    costDug = min(P.costDug_a*exp(P.costDug_b*lift),10*P.dDugInt); %limit the cost to keep absurdly high costs from blowing up the problem
    costBore = P.electricity*lift;
    newShr = invest + shrBore;
    
%     icPts = numel(P.investCosts);
% 
%     %find pseudo parcels that switched;
%     investCosts = repmat(P.investCosts,ns,1,3);
%     shrPts = repmat(P.shrPoints,ns,1,3);
%     probs = repmat(P.icProbs,ns,1,3);
%     newShrMat = repmat(newShr,1,icPts);
%     
%     deltaInvest = 1e-3;
%     oldShrMat = repmat(shrBore,1,icPts,3);
%     newShrMats(:,:,1) = newShrMat-deltaInvest;
%     newShrMats(:,:,2) = newShrMat;
%     newShrMats(:,:,3) = newShrMat+deltaInvest;
%     switched = (shrPts>oldShrMat).*(shrPts<=newShrMats);
%     investCostMat = squeeze(sum(investCosts.*switched.*probs,2));
%     
%     dinvestCost_dinvest = (investCostMat(:,3)-investCostMat(:,1))/(2*deltaInvest);
%     ddinvestCost_ddinvest = (investCostMat(:,3) + investCostMat(:,1) - 2*investCostMat(:,2))./(deltaInvest^2);
%     
% % %approximate the derivative by finding the closest point and calculating the cost
%     dShr = max(0,newShrMat-shrPts);
%     [a,ptSub] = min(dShr,[],2);
%     nextPtSub = min(icPts,ptSub+1);
%     ptInd = sub2ind([ns icPts],(1:ns)',ptSub);
%     nextPt = sub2ind([ns icPts],(1:ns)',nextPtSub);
%     dinvestCost_dinvest = investCosts(ptInd).*probs(ptInd)/(P.shrPoints(2)-P.shrPoints(1));
%     ddinvestCost_ddinvest = (investCosts(nextPt)-investCosts(ptInd)).*probs(ptInd)/((P.shrPoints(2)-P.shrPoints(1))^2);
%     
% %     %investment cost by parcel is distributed normally with mean investMean
%     %and variance investVar
%     lowCost = min(max(norminv(shrBore,P.investCostMean,P.investCostSD),-1/eps),1/eps); %Functionally truncate the sample
%     highCost = min(max(norminv(newShr,P.investCostMean,P.investCostSD),-1/eps),1/eps);
%     deltaCost = highCost-lowCost;
% %   

lowCost = norminv(shrBore,P.investCostMean,P.investCostSD);
highCost = norminv(newShr,P.investCostMean,P.investCostSD);
%compute the mean as the mean of the normal distribution truncated to the
%range (lowCost, highCost) times the probability a parcel lies in that
%range (e.g. F(highCost) - F(lowCost))
% 
%     approxShare = repmat([0:.0001:1],ns,1);
%     costArray = repmat(lowCost,1,size(approxShare,2)) + repmat(deltaCost,1,size(approxShare,2)).*approxShare;
%     probArray = normpdf(costArray,P.investCostMean,P.investCostSD);
%     intervalIntegral = sum(costArray.*probArray,2);
    
%     intervalIntegral = zeros(numel(lowCost),1);
%     for ii=1:numel(lowCost);
%         intervalIntegral(ii) = integral(@(x) x.*normpdf(x,P.investCostMean,P.investCostSD),lowCost(ii),highCost(ii));
%     end
% %     
    
% dinvestCost_dinvest = intervalIntegral + highCost.*invest;
% ddinvestCost_ddinvest = 1./normpdf(highCost,P.investCostMean,P.investCostSD)+2*highCost;
 
investCost = P.investCostMean*(newShr-shrBore) + P.investCostSD*(normpdf(lowCost,P.investCostMean,P.investCostSD)-normpdf(highCost,P.investCostMean,P.investCostSD)); %expected cost conditional on being in interval
dinvestCost_dinvest = P.investCostMean*(1-1/P.investCostSD) + highCost/P.investCostSD;
ddinvestCost_ddinvest = 1./(P.investCostSD*normpdf(highCost,P.investCostMean,P.investCostSD));
%     if any(invest); keyboard; end
    
%      keyboard
	nbDug = P.dDugInt*gwDug - P.dDugSlope/2*gwDug.^2 - costDug.*gwDug;
    nbBore = P.dBoreInt*gwBore - P.dBoreSlope/2*gwBore.^2 - costBore.*gwBore;
    
    nb.all = nbDug.*(1-shrBore) + nbBore.*shrBore - investCost;
    nb.dug = nbDug;
    nb.bore = nbBore;
    
%     if  max(abs(nb.all))>=100;
%         keyboard;
%     end
    
   	dnb.dgwDug = (P.dDugInt-P.dDugSlope.*gwDug-costDug).*(1-shrBore);
    dnb.dgwBore = (P.dBoreInt-P.dBoreSlope.*gwBore-costBore).*shrBore;
    dnb.di = -dinvestCost_dinvest;
    ddnb.ddgwDug = -P.dDugSlope*ones(size(gwLevel)).*(1-shrBore);
    ddnb.ddgwBore = -P.dBoreSlope*ones(size(gwLevel)).*shrBore;
    ddnb.dii = -ddinvestCost_ddinvest;
    
    %no cross partials
    ddnb.dgwDugdgwBore = zeros(size(gwLevel));
    ddnb.dgwDugdi = zeros(size(gwLevel));
    ddnb.dgwBoredi = zeros(size(gwLevel));
 