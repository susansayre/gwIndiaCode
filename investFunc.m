function out = investFunc(shares,investActions,P,flag)

%shares is ns x P.numTech and is the share of farms in each technology
%investActions is ns x P.numTech*(P.numTech-1)
%first P.numTech-1 are the share of farms currently in tech1 ending in tech 1, share of farms currently in tech 1 and
%not ending in farm 1 ending in farm 2 and so on
%second P.numTech-1 are the share of farms currently in tech2 ending tech 1, share of farms currently in tech 2 and not
%ending in farm 1 ending in farm 2 and so on

investActions = min(max(investActions,0),1);
ns = size(shares,1);

F1 = reshape(investActions,[ns P.numTech-1 P.numTech]);        
F1(:,P.numTech,:) = 1;

sharesArray = permute(repmat(shares,[1 1 P.numTech]),[1 3 2]); %ns x possibleStates
newSharesArray(1:ns,1,:) = sharesArray(:,1,:).*F1(:,1,:);
for ii=2:P.numTech
    newSharesArray(1:ns,ii,:) = prod(1-F1(:,1:ii-1,:),2).*F1(:,ii,:).*sharesArray(:,ii,:);
end

switch flag
    case 'tech'
        newShares = sum(newSharesArray,3);
        newSharesRaw = zeros(ns,P.numTech-1);
        newSharesRaw(:,1) =newShares(:,1);
        for ii=2:P.numTech-1
            myInds = intersect(find(newShares(:,ii)>0),find(sum(newShares(:,1:ii-1),2)<1));           
            newSharesRaw(myInds,ii) = newShares(myInds,ii)./(1-sum(newShares(myInds,1:ii-1),2));           
        end
        out = newSharesRaw;
        
    case 'cost'
        investAmtsHat = permute(newSharesArray,[3 2 1]);
        investBaseCosts = investAmtsHat.*repmat(P.techCostMat,[1 1 ns]);
        nsBigMat =permute(repmat((1:ns)',[1 P.numTech P.numTech]),[2 3 1]);
        orderedInvestAmts = investAmtsHat(sub2ind([P.numTech P.numTech ns],repmat((1:P.numTech)',1,P.numTech,ns),repmat(P.techCostOrder,1,1,ns),nsBigMat));
        orderedInvestCosts = investBaseCosts(sub2ind([P.numTech P.numTech ns],repmat((1:P.numTech)',1,P.numTech,ns),repmat(P.techCostOrder,1,1,ns),nsBigMat));

        cumShareEnd = permute(cumsum(sharesArray,3,'reverse')-sharesArray,[3 2 1]) + cumsum(orderedInvestAmts,2);
        cumShareStart = permute(cumsum(sharesArray,3,'reverse')-sharesArray,[3 2 1]);
        cumShareStart(:,2:P.numTech,:) = cumShareEnd(:,1:P.numTech-1,:);

        costStart = P.maxICM*cumShareStart;
        costEnd = P.maxICM*cumShareEnd;

        investCost = orderedInvestCosts.*(1+(costStart.*orderedInvestAmts+0.5*(costEnd-costStart).*orderedInvestAmts));  
 
        out =squeeze(sum(sum(investCost))) ;
end

if any(any(isnan(out)))||any(any(isinf(out)))
    if P.interactive
        keyboard
    else
        disp('Returning nans or infs in investFunc')
    end
end
   