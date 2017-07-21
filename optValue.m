function [value,valueDetail,farmDetail] = optValue(s,x,e,P,toWho)

if ~exist('toWho','var')
    toWho = 'society';
end

ns = size(s,1);

%modify parameters according to whose benefits we're computing
switch toWho
    case 'society'
        energyCosts = P.eCosts;
        fixedEcosts = repmat(P.FEcostSocShr.*P.fixedEcosts,ns,1);
    case 'indiv'
        energyCosts = P.eCosts*P.eCostShr;
        fixedEcosts = repmat(P.fixedEcosts,ns,1);
end

shares = rawShr2Share(s(:,P.shareInds));

lifts = P.landHeight - s(:,P.ind.level); %ns x 1
water = x(:,P.ind.water);

if any(e)
    error('This function can''t handle shocks')
end
 
costs = lifts*energyCosts; %ns x ntech

Intercepts = repmat(P.idInts,ns,1); %ns x ntech 
Slopes = repmat(P.idSlopes,ns,1); %ns x ntech

farmRev = Intercepts.*water - 0.5*Slopes.*(water.^2);
nb = farmRev - costs.*water - fixedEcosts;
investCost = investFunc(shares,x(:,P.ind.invest),P,'cost'); %ns x 1

value = sum(nb.*shares,2) - investCost;

if nargout>1
    valueDetail(:,1) = sum(farmRev.*shares,2);
    valueDetail(:,2) = sum(costs.*water.*shares,2);
    valueDetail(:,3) = sum(fixedEcosts.*shares,2);
    valueDetail(:,4) = investCost;
end

if nargout>2
    farmDetail.rev = farmRev;
    farmDetail.waterCost = costs.*water;
    farmDetail.tariff = fixedEcosts;
    farmDetail.profit = farmRev - costs.*water - fixedEcosts;
end

if any(isnan(value))
	if P.interactive
	    keyboard
	else
		disp('Returning nans in optValue')
	end
end

