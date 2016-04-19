function [out1,out2,out3] = optStoppingFunc(flag,s,x,e,P)
	
%states are water levels, percentile rank among farms, and type
%actions are invest (x=1) or don't invest (x=0)

    ns = size(s,1);
    ds = size(s,2);
    dx = size(x,2);

    switch flag
	
		case 'b'
            out1 = 0; %lower bounds on all variables are 0;
            out2 = 1;
            
        case 'f'
            %return net benefits
            levels = s(:,1);
            percCheaper = s(:,2); %state variable runs from 0 to 1 and identifies farm costs based on where in the interval the farm ranks
            costs = min(max(norminv(percCheaper,P.investCostMean,P.investCostSD),-1/eps),1/eps);
            type = s(:,3); %0=dug, 1=bore
            
            lift = P.landHeight - levels;
            costDug = P.costDug_a*exp(P.costDug_b*lift);
            costBore = P.electricity*lift;

            gwDug = (levels>0).*max(0,(P.dDugInt-costDug)/P.dDugSlope);
            gwBore = (levels>0).*max(0,(P.dBoreInt-costBore)/P.dBoreSlope);
    
            nbDug = P.dDugInt*gwDug - P.dDugSlope/2*gwDug.^2 - costDug.*gwDug;
            nbBore = P.dBoreInt*gwBore - P.dBoreSlope/2*gwBore.^2 - costBore.*gwBore;
            
            out1 = repmat((1-type).*nbDug + type.*nbBore,1,dx) - x.*repmat(costs,1,dx);
            
 		case 'g'
            %state transition
            out1(:,1) = max(P.bottom,s(:,1)+P.levelTrend);
            out1(:,2) = s(:,2);
            out1(:,3) = min(1,s(:,3)+x);
			
	end
		