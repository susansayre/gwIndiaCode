function [out1,out2,out3] = cpFunc(flag,s,x,e,P)
	shrBore = s(:,1);
	investmentAmts = x(:,P.investInd);
	gwDug = max(0,x(:,P.gwDugInd));
    gwBore = max(0,x(:,P.gwBoreInd));
	yrsLeft = s(:,2);

    if ~isfield(P,'explicitPath')
        P.explicitPath = 0;
    end
    if P.explicitPath
        maxYrsLeft = length(P.levelPath);
        if any(yrsLeft>maxYrsLeft)
            keyboard
        end
        gwLevels = 0*yrsLeft;
        evenYears = round(yrsLeft);
        gwLevels(find(evenYears)) = P.levelPath(maxYrsLeft+1-evenYears(find(evenYears)));
    else
        gwLevels = min(P.landHeight,max(0,polyval(P.levelParams,yrsLeft)));
    end

    if any(e)
        error(['This function can''t handle shocks'])
    end
    
    ns = size(s,1);
    ds = size(s,2);
    dx = size(x,2);

    switch flag
	
		case 'b'
            out1 = zeros(ns,dx); %lower bounds on all variables are 0;
            out2 = zeros(ns,dx);
            out2(:,P.gwDugInd) = P.dDugInt/P.dDugSlope; %we know the benefit is negative at higher amts
            out2(:,P.gwBoreInd) = P.dBoreInt/P.dBoreSlope; %we know the benefit is negative at higher amts
            out2(:,P.investInd) = min(.1,(1-s(:,P.sbInd)));  %this implies converting all of the additional parcels

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
			newYrsLeft = max(0,yrsLeft - 1);
			out1 = [shrBore+investmentAmts newYrsLeft];
            out2 = zeros(ns,ds,dx);
			out2(:,1,P.investInd) = ones(ns,1);
			out3 = zeros(ns,ds,dx,dx);
	end
		