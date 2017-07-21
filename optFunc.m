function [out1,out2,out3] = optFunc(flag,s,x,e,P)
  
	switch flag
	
		case 'b'
            ns = size(s,1);
            lifts = P.landHeight - s(:,P.ind.level); %ns x 1
            %return upper and lower bounds on the control variables
            out1 = zeros(ns,P.numTech^2); %lower bounds on all variables are 0;
            maxWater = waterLimit(lifts,P);
            
            shares = rawShr2Share(s(:,P.shareInds));
            maxWater(shares==0) = 0;
            
            out2(:,P.ind.water) = maxWater;
            if P.maxInvest==0
                investLimits = eye(P.numTech);
                maxInvest = repmat(reshape(investLimits(:,1:P.numTech-1),1,P.numTech*(P.numTech-1)),ns,1);
            elseif P.maxInvest==1
                maxInvest = P.maxInvest*ones(ns,P.numTech*(P.numTech-1));
            else
                keyboard
            end
            
            for ii=1:P.numTech
                maxInvest(shares(:,ii)==0,(ii-1)*(P.numTech-1)+1:ii*(P.numTech-1)) = 0; %if there aren't any farms with a technology don't allow any movement from that technology
            end
            out2(:,P.ind.invest) = maxInvest;


   		case 'f'
            out1 = optValue(s,x,e,P);
            if nargout > 1
                out2 = fd(s,x,e,P,'optValue');
            end
            if nargout > 2
                out3 = fd(s,x,e,P,'optValue','fd');
            end
            
        case 'g'
            %return updated state
            out1 = optState(s,x,e,P);
            if nargout > 1
                out2 = fd(s,x,e,P,'optState');
            end
            if nargout > 2
                out3 = fd(s,x,e,P,'optState','fd');
            end
	end
	
    if any(isnan(out1))
		if P.interactive
			keyboard
		else
			disp('Returning nans in optFunc')
		end
    end