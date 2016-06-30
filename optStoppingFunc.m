function [out1,out2,out3] = optStoppingFunc(flag,s,x,e,P,liftParams)
	
%states are water levels, percentile rank among farms, and type
%actions are invest (x=1) or don't invest (x=0)

    ns = size(s,1);
    ds = size(s,2);
    dx = size(x,2);

    switch flag
	
		case 'b'
            out1 = 0*x; %lower bounds on all variables are 0;
            out2 = 1;
            
        case 'f'
            %return net benefits
            levels = s(:,1);
            percCheaper = s(:,2); %state variable runs from 0 to 1 and identifies farm costs based on where in the interval the farm ranks
            costs = norminv(percCheaper,P.investCostMean,P.investCostSD);
            type = s(:,3); %0=dug, 1=bore
            
            lift = P.landHeight - levels;
            costDug = P.costDug_a*exp(P.costDug_b*lift);
            costBore = P.electricity*lift; %+ exp(P.costBore_a+P.costBore_b*levels);

            gwDug = (levels>0).*max(0,(P.idDugInt-costDug)/P.idDugSlope);
            gwBore = (levels>0).*max(0,(P.idBoreInt-costBore)/P.idBoreSlope);
    
            nbDug = P.idDugInt*gwDug - P.idDugSlope/2*gwDug.^2 - costDug.*gwDug;
            nbBore = P.idBoreInt*gwBore - P.idBoreSlope/2*gwBore.^2 - costBore.*gwBore;
            
            out1 = min(max(repmat((1-type).*nbDug + type.*nbBore,1,dx) - x.*repmat(costs+P.convertTax,1,dx),-1/eps),1/eps);
            
 		case 'g'
            %state transition
            lift = max(eps,P.landHeight - s(:,1));
            estTnext = liftParams.t0+1 - 1/(liftParams.growthRate*liftParams.nu)*(log(max(1+eps,liftParams.ss./lift).^liftParams.nu-1));
            liftNext = liftParams.ss./((1+exp(-liftParams.growthRate*liftParams.nu*(estTnext-liftParams.t0))).^(1/liftParams.nu));
            out1(:,1) = max(P.bottom,P.landHeight-liftNext);
            out1(:,2) = s(:,2);
            out1(:,3) = min(1,s(:,3)+x);
            if ~isreal(out1), keyboard, end
    end

    if any(isnan(out1));
        keyboard
    end