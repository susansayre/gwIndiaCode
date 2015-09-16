function [out1,out2,out3] = optFunc(flag,s,x,e,P)
	gwLevels = s(:,2);
	wellCapitalLevels = s(:,1);
	investmentAmts = x(:,P.investInd);
	gwExtractionAmts = x(:,P.gwInd);
	
    if any(e)
        error(['This function can''t handle shocks'])
    end
    
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
			[g,dg,dgg] = updateLevels(gwLevels,gwExtractionAmts,P);
			newWellCap = (1-P.wellCapDecayRate)*wellCapitalLevels + investmentAmts;
			out1 = [newWellCap g];
			out2(:,1,P.investInd) = ones(ns,1); %extraction doesn't affect capital
            out2(:,2,P.gwInd) = dg; %investment doesn't affect level
			out3(:,2,P.gwInd,P.gwInd) = dgg; %likely zero but set by updateLevels function
	end
		