function [out1,out2,out3] = optFunc(flag,s,x,e,P)
	gwLevels = s(:,P.levelInd);
	wellCapitalLevels = s(:,P.wellInd);
	investmentAmts = x(:,P.investInd);
	gwExtractionAmts = x(:,P.gwInd);
	
    if any(e)
        error(['This function can''t handle shocks'])
    end
    
    ns = size(s,1);
    ds = size(s,2);
    dx = size(x,2);
    
	switch flag
	
		case 'b'
            %return upper and lower bounds on the control variables
            out1 = repmat([0 0],ns,1);
            maxGW = P.gwdIntercept/P.gwdSlope;
            lift = P.landHeight - gwLevels;
            gwUBDug = (lift<=P.maxDepthNoCap)*P.noCapLimitShare*maxGW;
            out2(:,P.gwInd) = wellCapitalLevels*maxGW+(1-wellCapitalLevels).*gwUBDug; %above this, the net benefit of gw is negative
            out2(:,P.investInd) = P.maxInvest; 

   		case 'f'
            %return net benefits
			[out1, dnb, ddnb] = netBen(gwExtractionAmts,investmentAmts,gwLevels,wellCapitalLevels,P);
            out2(:,P.investInd) = dnb.di;
            out2(:,P.gwInd) = dnb.dgw;
            out3(:,P.investInd,P.investInd)=ddnb.dii;
            out3(:,P.investInd,P.gwInd) = ddnb.dgwdi;
            out3(:,P.gwInd,P.investInd) = ddnb.dgwdi;
            out3(:,P.gwInd,P.gwInd) = ddnb.dgwgw;
        case 'g'
            %return updated states
			[g,dg,dgg] = updateLevels(gwLevels,gwExtractionAmts,P);
			newWellCap = (1-P.wellCapDecayRate)*wellCapitalLevels + investmentAmts;
			out1(:,P.wellInd) = newWellCap; 
            out1(:,P.levelInd) = g;
            
            %return derivatives of next period states with respect to actions
            out2(:,P.wellInd,P.investInd) = ones(ns,1); 
            out2(:,P.wellInd,P.gwInd) = zeros(ns,1); %extraction doesn't affect capital
            out2(:,P.levelInd,P.gwInd) = dg;
            out2(:,P.levelInd,P.investInd) = zeros(ns,1); %investment doesn't affect levels directly
			
            %return second derivatives of next period states with respect to actions
            out3 = zeros(ns,ds,dx,dx); %initialize to zero and add only as needed
            out3(:,P.levelInd,P.gwInd,P.gwInd) = dgg; %likely zero but set by updateLevels function
	end
		