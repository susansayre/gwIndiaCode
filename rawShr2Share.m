function shares = rawShr2Share(rawShares)

ns = size(rawShares,1);
nT = size(rawShares,2) + 1;
shares = zeros(ns,nT);
shares(:,1) = rawShares(:,1);
for ii=2:nT-1
    shares(:,ii) = (1-sum(shares(:,1:ii-1),2)).*rawShares(:,ii);
end
shares(:,nT) = 1-sum(shares(:,1:nT-1),2);
