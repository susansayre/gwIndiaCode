function [nb,dnb,ddnb] = netBen(gw,invest,gwLevel,capital,P)
	lift = P.landHeight - gwLevel;
	alpha = 0 + P.electricity*capital; 
    beta = P.gwdIntercept/(P.maxDepthNoCap^2)*(1-capital); %when capital is zero, cost per unit is greater than benefit; when capital is one, there is no quadratic term
	beta = 0; %try alternate specification
    costPerUnit = alpha.*lift + beta.*(lift.^2); %per volume cost of gw as a function of depth and well capital
	investCost = P.fixedInvestCost*(invest>0) + P.varInvestCost*invest; %cost of increasing well capital as a function of investment
	nb = P.gwdIntercept*gw - P.gwdSlope/2*gw.^2 - costPerUnit.*gw - investCost;
    if max(nb)>10e5
        keyboard
    end   
	
	dnb.dgw = P.gwdIntercept - P.gwdSlope*gw - costPerUnit;
    dnb.di = -P.varInvestCost*ones(size(gwLevel));
    ddnb.dgwgw = -P.gwdSlope*ones(size(gwLevel));
    ddnb.dgwdi = zeros(size(gwLevel));
    ddnb.dii = zeros(size(gwLevel));
    
 