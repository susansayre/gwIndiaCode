function aeOutput = aeSolve2(P,modelOpts)
%solve the adaptive expectations problem as an optimal switching problem

% statePath(1,P.levelInd) = P.h0; 
% statePath(1,P.sbInd) = P.shrBore0;
% trendPath = P.levelTrend;

t=1;
% change = 1;
P.iTol = 1e-3;

model.func = 'optStoppingFunc';
model.discount = P.discount;
model.actions =  [0;1];
model.discretestates = 3;

model.e = 0;
model.w = 1;


levelPath = P.h0;
shrBore = P.shrBore0;
% trend = P.levelTrend;
valChange = 1;

iter = 0;
while valChange>P.iTol
    iter = iter+1;
    %find maxCost of farms that have already adopted.
    lowCost = norminv(shrBore(t)*P.inTruncProb+P.probBelow,P.investCostMean,P.investCostSD);
    %solve optimal stopping problem for this trend
    %convert state variable to be naturally bounded on [0,1].
    fspace = fundefn('lin',[modelOpts.heightNodes modelOpts.capNodes],[P.bottom shrBore(t)],[P.landHeight 1],[],[0;1]);
    scoord = funnode(fspace);
    snodes = gridmake(scoord);
    model.params = {P};
    optset('dpsolve','showiters',0)
    figure()
    [c,s,v,x] = dpsolve(model,fspace,snodes);
    close
    %simulate forward to see who invests this period
    nyrs = 1;
    
    simulStates = [repmat(levelPath(t),numel(s{2}),1) s{2} 0*s{2}];
    [spath,xpath] = dpsimul(model,simulStates,nyrs,s,x);
    %find maxCost of farms that have adopted by end of period
    shrBore(t+1) = max(shrBore(t),max(spath(:,3,2).*spath(:,2,2)));
    invest = shrBore(t+1) - shrBore(t);
    highCost = norminv(shrBore(t+1)*P.inTruncProb + P.probBelow,P.investCostMean,P.investCostSD);
    if highCost<=lowCost
        investCost = 0;
    else
        investCost = P.investCostMean*invest+P.investCostSD^2/P.inTruncProb*(normpdf(lowCost,P.investCostMean,P.investCostSD)-normpdf(highCost,P.investCostMean,P.investCostSD));
    end
    
    lift = P.landHeight - levelPath(t);
    costDug = P.costDug_a*exp(P.costDug_b*lift);
    costBore = P.electricity*lift;

    gwDug = (levelPath(t)>0).*max(0,(P.idDugInt-costDug)/P.idDugSlope);
    gwBore = (levelPath(t)>0).*max(0,(P.idBoreInt-costBore)/P.idBoreSlope);

    nbDug = P.idDugInt*gwDug - P.idDugSlope/2*gwDug.^2 - costDug.*gwDug;
    nbBore = P.idBoreInt*gwBore - P.idBoreSlope/2*gwBore.^2 - costBore.*gwBore;
    
    gwUse = (1-shrBore(t))*gwDug + shrBore(t)*gwBore;
    
    levelPath(t+1) = updateLevels(levelPath(t),gwUse,P);
    P.levelTrend = levelPath(t+1)-levelPath(t);
    
    val(t,:) = [nbDug nbBore (1-shrBore(t))*nbDug + shrBore(t)*nbBore - investCost invest*P.convertTax];
    xPath(t,:) = [shrBore(t+1)-shrBore(t) gwDug gwBore];

    valChange = abs(val(t,3)*P.discount^(t-1));
    fprintf ('%4i %10.1e\n',iter,valChange)
    t = t+1;

end

aeOutput.statePath = [shrBore' levelPath'];
aeOutput.controlPath = xPath;
aeOutput.valPath = val;
aeOutput.aeVal = (P.discount.^(0:length(val)-1))*val(:,3);
   
	
