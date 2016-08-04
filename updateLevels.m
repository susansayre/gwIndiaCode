function [newLevel,dgw,ddgw] = updateLevels(oldLevel,gw,P)
	netRecharge = P.inflow - (1-P.returnFlow)*gw;
	newLevel = oldLevel + 1/P.AS*netRecharge;
    atTop = (newLevel>P.landHeight-.5);
    newLevel(find(atTop)) = P.landHeight-P.minLift;
	dgw = -1/P.AS*(1-P.returnFlow)*ones(size(gw));
    dgw(find(atTop)) = 0;
	ddgw = zeros(size(gw));