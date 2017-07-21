function maxWater = waterLimit(lifts,P,ind)
%takes a vector of pumping lifts and returns the maximum water extractable from each technology

ns = size(lifts,1);

fracTowardBottom = 1-min(1,max(repmat(P.maxDepths,ns,1) - repmat(lifts,1,P.numTech),0)./repmat(P.maxDepths - P.maxLiftMaxCap,ns,1));
%share = 1./(1+exp(P.logisticRate*(fracTowardBottom-0.5)));
share = 1-fracTowardBottom;
% share(fracTowardBottom==1) =0;
% share(fracTowardBottom==0) = 1;

maxWater = repmat(P.maxWater*P.maxCapacity,ns,1).*share;

if exist('ind','var')
    maxWater = maxWater(ns*(ind-1)+(1:ns)');
end