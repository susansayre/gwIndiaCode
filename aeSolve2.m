function aeOutput = aeSolve2(P,modelOpts)
%solve the adaptive expectations problem as an optimal switching problem
numYrs = 50;
statePath(1,P.levelInd) = P.h0; 
statePath(1,P.sbInd) = P.shrBore0;
trendPath = P.levelTrend;

t=1;
change = 1;
P.iTol = 1e-5;

model.func = 'optStoppingFunc';
model.discount = P.discount;
model.actions =  [0;1];
model.discretestates = 3;

model.e = 0;
model.w = 1;


levelPath = P.h0;
shrBore = P.shrBore0;
trend = P.levelTrend;
valChange = 1;

while valChange>P.iTol
    
    %find maxCost of farms that have already adopted.

    lowCost = min(max(norminv(shrBore(t),P.investCostMean,P.investCostSD),-1/eps),1/eps);
    %solve optimal stopping problem for this trend
    fspace = fundefn('spli',[10 10],[P.bottom lowCost],[P.landHeight P.investCostMean+3*P.investCostSD],[],[0;1]);
    scoord = funnode(fspace);
    snodes = gridmake(scoord);
    model.params = {P};
    [c,s,v,x] = dpsolve(model,fspace,snodes);

    %simulate forward to see who invests this period
    nyrs = 1;
    
    simulStates = [repmat(levelPath(t),numel(s{2}),1) s{2} 0*s{2}];
    [spath,xpath] = dpsimul(model,simulStates,1,s,x);
    %find maxCost of farms that have adopted by end of period
    highCost = max(spath(:,3,2).*spath(:,2,2));
    if highCost<=lowCost
        shrBore(t+1) = shrBore(t);
        investCost = 0;
    else
        shrBore(t+1) = normcdf(highCost,P.investCostMean,P.investCostSD);
        investCost = integral(@(x) x.*normpdf(x,P.investCostMean,P.investCostSD),lowCost,highCost)*(shrBore(t+1)-shrBore(t));
    end
    
    lift = P.landHeight - levelPath(t);
    costDug = P.costDug_a*exp(P.costDug_b*lift);
    costBore = P.electricity*lift;

    gwDug = (levelPath(t)>0).*max(0,(P.dDugInt-costDug)/P.dDugSlope);
    gwBore = (levelPath(t)>0).*max(0,(P.dBoreInt-costBore)/P.dBoreSlope);

    nbDug = P.dDugInt*gwDug - P.dDugSlope/2*gwDug.^2 - costDug.*gwDug;
    nbBore = P.dBoreInt*gwBore - P.dBoreSlope/2*gwBore.^2 - costBore.*gwBore;
    
    gwUse = (1-shrBore(t))*gwDug + shrBore(t)*gwBore;
    
    levelPath(t+1) = updateLevels(levelPath(t),gwUse,P);
    P.levelTrend = levelPath(t+1)-levelPath(t);
    
    val(t,:) = [nbDug nbBore (1-shrBore(t))*nbDug + shrBore(t)*nbBore - investCost];
    xPath(t,:) = [shrBore(t+1)-shrBore(t) gwDug gwBore];

    valChange = val(t,3)*P.discount^(t-1);
    t = t+1;

end

aeOutput.statePath = [shrBore' levelPath'];
aeOutput.controlPath = xPath;
aeOutput.valPath = val;
aeOutput.aeVal = P.discount.^(0:length(val)-1)*val(:,3);
   
	
