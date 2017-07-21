function sharesRaw = share2RawShr(shares)

[ns,nT] = size(shares);
sharesRaw = zeros(ns,nT-1);
sharesRaw(:,1) = shares(:,1);
for ii=2:nT-1
    myInds = find(shares(:,ii)>0);
    sharesRaw(myInds,ii) = shares(myInds,ii)./(1-sum(shares(myInds,1:ii-1),2));           
end