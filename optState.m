function newState = optState(s,x,e,P)

ns = size(s,1);

rawShares = s(:,P.shareInds);
shares = zeros(ns,P.numTech);
shares(:,1) = rawShares(:,1);
for ii=2:P.numTech-1
    shares(:,ii) = (1-sum(shares(:,1:ii-1),2)).*rawShares(:,ii);
end
shares(:,P.numTech) = 1-sum(shares(:,1:P.numTech-1),2);

water = x(:,P.ind.water);

if any(e)
    error('This function can''t handle shocks')
end
 
gwExtractionAmts = sum(water.*shares,2);
g = updateLevels(s(:,P.ind.level),gwExtractionAmts,P);

newState(:,P.ind.level) = g;
newState(:,P.shareInds) = investFunc(shares,x(:,P.ind.invest),P,'tech');

if any(any(isnan(newState)))
    disp('Returning nans in optState')
end