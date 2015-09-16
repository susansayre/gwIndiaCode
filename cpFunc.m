function [out1,out2,out3] = cpFunc(flag,s,x,e,P)
	wellCapitalLevels = s(:,1);
	investmentAmts = x(:,P.investInd);
	gwExtractionAmts = x(:,P.gwInd);
	yrsLeft = s(:,2);
	gwLevels = polyval(P.levelParams,yrsLeft);
	
    ns = size(s,1);
	switch flag
	
		case 'b'
	        out1 = repmat([0 0],ns,1);
            out2 = repmat([P.gwdIntercept/P.gwdSlope P.maxInvest],ns,1);
		case 'f'
			[out1, dnb, ddnb] = netBen(gwExtractionAmts,investmentAmts,gwLevels,wellCapitalLevels,P);
            out2(:,P.investInd) = dnb.di;
            out2(:,P.gwInd) = dnb.dgw;
            out3(:,P.investInd,P.investInd)=ddnb.dii;
            out3(:,P.investInd,P.gwInd) = ddnb.dgwdi;
            out3(:,P.gwInd,P.investInd) = ddnb.dgwdi;
            out3(:,P.gwInd,P.gwInd) = ddnb.dgwgw;
 		case 'g'
			newWellCap = (1-P.wellCapDecayRate)*wellCapitalLevels + investmentAmts;
			newYrsLeft = max(0,yrsLeft - 1);
			out1 = [newWellCap newYrsLeft];
            out2 = zeros(ns,2,2);
			out2(:,1,P.investInd) = ones(ns,1);
			out3 = zeros(ns,size(s,2),size(x,2),size(x,2));
	end
		