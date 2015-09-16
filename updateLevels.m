function [newLevel,dgw,ddgw] = updateLevels(oldLevel,gw,P)
	netRecharge = P.inflow - (1-P.returnFlow)*gw;
	newLevel = oldLevel + 1/P.AS*netRecharge;
	dgw = -1/P.AS*(1-P.returnFlow)*ones(size(gw));
	ddgw = zeros(size(gw));