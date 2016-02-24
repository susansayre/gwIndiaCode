function [nb,dnb,ddnb] = netBen(gwDug,gwBore,invest,gwLevel,shrBore,P)
    %this function returns the aggregate net benefit as a function of gw
    %extraction (which now needs to be a two dimensional vector),
    %investment in regime switching, gwLevels, and the current share of
    %plots that are irrigated by bore wells;
    
    ns = length(shrBore);
    
    gwLevel = max(P.bottom,min(P.landHeight,gwLevel));
	lift = P.landHeight - gwLevel;
    
    if any(lift>1)||any(lift<0)
        keyboard
    end
    
    %compute cost per unit for dug wells
    costDug = P.costDug_a*exp(P.costDug_b*lift);
    costBore = P.electricity*lift;
    newShr = invest + shrBore;
     
    %investment cost by parcel is distributed normally with mean investMean
    %and variance investVar
    lowCost = min(max(norminv(shrBore,P.investCostMean,P.investCostSD),-1/eps),1/eps);
    highCost = min(max(norminv(newShr,P.investCostMean,P.investCostSD),-1/eps),1/eps);
    deltaCost = highCost-lowCost;
    
    approxShare = repmat([0:.1:1],ns,1);
    costArray = repmat(lowCost,1,size(approxShare,2)) + repmat(deltaCost,1,size(approxShare,2)).*approxShare;
    probArray = normpdf(costArray,P.investCostMean,P.investCostSD);
    intervalIntegral = sum(costArray.*probArray,2);
    
%     intervalIntegral = zeros(numel(lowCost),1);
%     for ii=1:numel(lowCost);
%         intervalIntegral(ii) = integral(@(x) x.*normpdf(x,P.investCostMean,P.investCostSD),lowCost(ii),highCost(ii));
%     end
    
    investCost = intervalIntegral.*invest; %expect cost conditional on being in interval
    
    dinvestCost_dinvest = intervalIntegral + highCost;
    ddinvestCost_ddinvest = 1./normpdf(highCost,P.investCostMean,P.investCostSD)+highCost;
    
	nbDug = P.dDugInt*gwDug - P.dDugSlope/2*gwDug.^2 - costDug.*gwDug;
    nbBore = P.dBoreInt*gwBore - P.dBoreSlope/2*gwBore.^2 - costBore.*gwBore;
    
    nb.all = nbDug.*(1-shrBore) + nbBore.*shrBore - investCost;
    nb.dug = nbDug;
    nb.bore = nbBore;
    
    if abs(max(nb.all))>=1e10;
        keyboard;
    end
    
   	dnb.dgwDug = P.dDugInt-P.dDugSlope.*gwDug-costDug;
    dnb.dgwBore = P.dBoreInt-P.dBoreSlope.*gwBore-costBore;
    dnb.di = -dinvestCost_dinvest;
    ddnb.ddgwDug = -P.dDugSlope*ones(size(gwLevel));
    ddnb.ddgwBore = -P.dBoreSlope*ones(size(gwLevel));
    ddnb.dii = -ddinvestCost_ddinvest;
    
    %no cross partials
    ddnb.dgwDugdgwBore = zeros(size(gwLevel));
    ddnb.dgwDugdi = zeros(size(gwLevel));
    ddnb.dgwBoredi = zeros(size(gwLevel));
 