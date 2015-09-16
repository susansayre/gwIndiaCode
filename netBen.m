function [nb,dnb,ddnb] = netBen(gw,invest,gwLevel,capital,P)
	lift = P.landHeight - gwLevel;
	alpha = 0 + P.electricity*capital; 
    beta = P.gwdIntercept/P.maxDepthNoCap^2*capital;
	costPerUnit = alpha.*lift + beta.*(lift.^2); %per volume cost of gw as a function of depth and well capital
	investCost = P.fixedInvestCost*(invest>0) + P.varInvestCost*invest; %cost of increasing well capital as a function of investment
	nb = P.gwdIntercept*gw - P.gwdSlope/2*gw.^2 - costPerUnit.*gw - investCost;
	
	dnb.dgw = P.gwdIntercept - P.gwdSlope*gw - costPerUnit;
    dnb.di = -P.varInvestCost*ones(size(gwLevel));
    ddnb.dgwgw = -P.gwdSlope*ones(size(gwLevel));
    ddnb.dgwdi = zeros(size(gwLevel));
    ddnb.dii = zeros(size(gwLevel));