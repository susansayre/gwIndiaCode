function [nb,dnb,ddnb] = netBen(gwDug,gwBore,invest,gwLevel,shrBore,P)
    %this function returns the aggregate net benefit as a function of gw
    %extraction (which now needs to be a two dimensional vector),
    %investment in regime switching, gwLevels, and the current share of
    %plots that are irrigated by bore wells;
    
    ns = length(shrBore);
    
   % gwLevel = max(P.bottom,min(P.landHeight,gwLevel));
	lift = P.landHeight - gwLevel;
    
    if any(lift>P.landHeight-P.bottom)||any(lift<0)
        keyboard
    end
    
    %compute cost per unit for dug wells
    costDug = P.electricityDug*lift;
    costBore = P.electricityBore*lift;
    newShr = invest + shrBore;
    

lowCost = norminv(shrBore*P.inTruncProb + P.probBelow,P.investCostMean,P.investCostSD);
highCost = norminv(min(newShr*P.inTruncProb  + P.probBelow,1),P.investCostMean,P.investCostSD);

aboveMaxInds = find(highCost>P.maxInvestCost); 
highCost(aboveMaxInds) = P.maxInvestCost;
investCost = P.investCostMean*invest + P.investCostSD^2*(normpdf(lowCost,P.investCostMean,P.investCostSD)-normpdf(highCost,P.investCostMean,P.investCostSD))/P.inTruncProb; %expected cost conditional on being in interval
dinvestCost_dinvest = highCost;
ddinvestCost_ddinvest = P.inTruncProb./normpdf(highCost,P.investCostMean,P.investCostSD);

    
%      keyboard
	nbDug = P.idDugInt*gwDug - P.idDugSlope/2*gwDug.^2 - costDug.*gwDug;
    nbBore = P.idBoreInt*gwBore - P.idBoreSlope/2*gwBore.^2 - costBore.*gwBore;
    
    nb.all = nbDug.*(1-shrBore) + nbBore.*shrBore - investCost - P.investCostPenalty*invest.^2;
    nb.dug = nbDug;
    nb.bore = nbBore;
    
    if  any(isnan(nb.all))
        keyboard;
    end
    
   	dnb.dgwDug = (P.idDugInt-P.idDugSlope.*gwDug-costDug).*(1-shrBore);
    dnb.dgwBore = (P.idBoreInt-P.idBoreSlope.*gwBore-costBore).*shrBore;
    dnb.di = -dinvestCost_dinvest -2*P.investCostPenalty*invest;
    ddnb.ddgwDug = -P.idDugSlope*ones(size(gwLevel)).*(1-shrBore);
    ddnb.ddgwBore = -P.idBoreSlope*ones(size(gwLevel)).*shrBore;
    ddnb.dii = -ddinvestCost_ddinvest - 2*P.investCostPenalty;
    
    %no cross partials
    ddnb.dgwDugdgwBore = zeros(size(gwLevel));
    ddnb.dgwDugdi = zeros(size(gwLevel));
    ddnb.dgwBoredi = zeros(size(gwLevel));
 