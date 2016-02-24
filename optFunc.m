function [out1,out2,out3] = optFunc(flag,s,x,e,P)
	gwLevels = s(:,P.levelInd);
	shrBore = s(:,P.sbInd);
	investmentAmts = x(:,P.investInd);
	gwDug = max(0,x(:,P.gwDugInd));
    gwBore = max(0,x(:,P.gwBoreInd));
	
    if any(e)
        error(['This function can''t handle shocks'])
    end
    
    ns = size(s,1);
    ds = size(s,2);
    dx = size(x,2);
    
	switch flag
	
		case 'b'
            %return upper and lower bounds on the control variables
            out1 = zeros(ns,dx); %lower bounds on all variables are 0;
            out2 = zeros(ns,dx);
            out2(:,P.gwDugInd) = P.dDugInt/P.dDugSlope; %we know the benefit is negative at higher amts
            out2(:,P.gwBoreInd) = P.dBoreInt/P.dBoreSlope; %we know the benefit is negative at higher amts
            out2(:,P.investInd) = (1-s(:,P.sbInd));  %this implies converting all of the additional parcels

   		case 'f'
            %return net benefits
			[b, dnb, ddnb] = netBen(gwDug,gwBore,investmentAmts,gwLevels,shrBore,P);
            out1 = b.all;
            out2(:,P.investInd) = dnb.di;
            out2(:,P.gwDugInd) = dnb.dgwDug;
            out2(:,P.gwBoreInd) = dnb.dgwBore;
            
            out3(:,P.investInd,P.investInd)=ddnb.dii;
            out3(:,P.investInd,P.gwDugInd) = ddnb.dgwDugdi;
            out3(:,P.investInd,P.gwBoreInd) = ddnb.dgwBoredi;
            
            out3(:,P.gwDugInd,P.gwDugInd) = ddnb.ddgwDug;
            out3(:,P.gwDugInd,P.gwBoreInd) = ddnb.dgwDugdgwBore;
            out3(:,P.gwDugInd,P.investInd) = ddnb.dgwDugdi;
            
            out3(:,P.gwBoreInd,P.gwBoreInd) = ddnb.ddgwBore;
            out3(:,P.gwBoreInd,P.gwDugInd) = ddnb.dgwDugdgwBore;
            out3(:,P.gwBoreInd,P.investInd) = ddnb.dgwBoredi;
            
        case 'g'
            %return updated states
            gwExtractionAmts = gwDug.*(1-shrBore) + gwBore.*shrBore;
			[g,dg,dgg] = updateLevels(gwLevels,gwExtractionAmts,P);
			out1(:,P.sbInd) = shrBore+investmentAmts;
            out1(:,P.levelInd) = g;
            
            %return derivatives of next period states with respect to actions
            out2(:,P.sbInd,P.investInd) = ones(ns,1); 
            out2(:,P.sbInd,P.gwDugInd)= zeros(ns,1); %extraction doesn't affect capital
            out2(:,P.sbInd,P.gwBoreInd)= zeros(ns,1); %extraction doesn't affect capital
            out2(:,P.levelInd,P.gwDugInd) = dg;
            out2(:,P.levelInd,P.gwBoreInd) = dg;
            out2(:,P.levelInd,P.investInd) = zeros(ns,1); %investment doesn't affect levels directly
			
            %return second derivatives of next period states with respect to actions
            out3 = zeros(ns,ds,dx,dx); %initialize to zero and add only as needed
            if any(dgg)
                error('have nto coded placement of second derivatives for gw levels')
            end
	end
		