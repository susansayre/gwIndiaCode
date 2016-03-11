function reOutput = reSolve(P,modelOpts)

t = 1;
model.func = 'optStoppingFunc';
model.discount = P.discount;
model.actions = [0;1]
model.discretestates = 3;
model.e = 0;
model.w = 1;

P.iTol = 1e-5;
P.trendTol = 1e-3;
shrBore = P.shrBore0;
valChange = 1;

levelPath = P.h0;
shrBore = P.shrBore0;

% optset('dpsolve','algorithm','funcit');
while valChange>P.iTol
    %find maxCost of farms that have already adopted.
    lowCost = min(max(norminv(shrBore(t),P.investCostMean,P.investCostSD),-1/eps),1/eps);
    %solve optimal stopping problem for this trend
    fspace = fundefn('spli',[10 10],[P.bottom lowCost],[P.landHeight P.investCostMean+3*P.investCostSD],[],[0;1]);
    scoord = funnode(fspace);
    snodes = gridmake(scoord);
    
    iter = 0;
    trendChange = 1;
    while trendChange>P.trendTol
        iter = iter + 1;
        model.params = {P};
        [c,s,v,x] = dpsolve(model,fspace,snodes);

        %simulate forward to see who invests under this trend  
        nyrs=1;

        simulStates = [repmat(levelPath(t),numel(s{2}),1) s{2} 0*s{2}];
        [spath,xpath] = dpsimul(model,simulStates,nyrs,s,x);
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

        newLevel = updateLevels(levelPath(t),gwUse,P);
        newTrend = newLevel -levelPath(t);
        
        trendChange = abs(P.levelTrend-newTrend)
        pause(3)
        P.levelTrend = newTrend;
    end
    
    levelPath(t+1) = newLevel;
    val(t,:) = [nbDug nbBore (1-shrBore(t))*nbDug + shrBore(t)*nbBore - investCost];
    xPath(t,:) = [shrBore(t+1)-shrBore(t) gwDug gwBore];

    valChange = val(t,3)*P.discount^(t-1);
    t = t+1;		
end

reOutput.statePath = [shrBore' levelPath'];
reOutput.controlPath = xPath;
reOutput.valPath = val;
reOutput.aeVal = P.discount.^(0:length(val)-1)*val(:,3);
