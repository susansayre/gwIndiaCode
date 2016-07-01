function reOutput = reSolve(P,modelOpts,initialPath)

t = 1;
model.func = 'optStoppingFunc';
model.discount = P.discount;
model.actions = [0;1];
model.discretestates = 3;
model.e = 0;
model.w = 1;

P.iTol = 1e-3;
P.trendTol = 1e-1;
shrBore = P.shrBore0;

assumedPath = initialPath;
pathChange = P.trendTol*10;
iter = 0;
cycleCount = 1;
moveFrac = .1;
% optset('dpsolve','algorithm','funcit');
while pathChange>P.trendTol
    t=1;
    levelPath = P.h0;
    shrBore = P.shrBore0;
    val = [];
    xPath = [];
    %solve optimal stopping problem for this trend
    optset('dpsolve','showiters',0)

    iter = iter + 1;

    fspace = fundefn('lin',[numel(assumedPath) modelOpts.capNodes*2],[1 shrBore(t)],[numel(assumedPath) 1],[],[0;1]);
    scoord = funnode(fspace);
    snodes = gridmake(scoord);

    model.params = {P 'path' assumedPath};
    [c,s,v,x] = dpsolve(model,fspace,snodes);
    close()
    %simulate forward to see who invests under this trend
    extraYrs = 40;
    assumedPath = [assumedPath; assumedPath(end)*ones(extraYrs,1)];
    nyrs = numel(assumedPath);
    simulStates = [ones(size(s{2})) s{2} 0*s{2}];
    [spath,xpath] = dpsimul(model,simulStates,numel(assumedPath),s,x);
    
    valChange = 10*P.iTol;
    while valChange>P.iTol
        %compute conversion shares and water use sequentially
 
        %find maxCost of farms that have already adopted.
        
        lift = P.landHeight - levelPath(t);
        costDug = P.costDug_a*exp(P.costDug_b*lift);
        costBore = P.electricity*lift;

        gwDug = (levelPath(t)>P.bottom).*max(0,(P.idDugInt-costDug)./P.idDugSlope);
        gwBore = (levelPath(t)>P.bottom).*max(0,(P.idBoreInt-costBore)./P.idBoreSlope);

        nbDug = P.idDugInt*gwDug - P.idDugSlope/2*gwDug.^2 - costDug.*gwDug;
        nbBore = P.idBoreInt*gwBore - P.idBoreSlope/2*gwBore.^2 - costBore.*gwBore;

        %spath(:,3,t+1) identifies who has converted by end of t+1
        %spath(:,2,t+1) gives the fraction of farms with costs <= a given farm
        %taking the max of the product of these gives the share cheaper than the last farm to convert, e.g. the
        %share converting
        lowCost = norminv(shrBore(t)*P.inTruncProb+P.probBelow,P.investCostMean,P.investCostSD);
        shrBore(t+1) = max(shrBore(t),max(spath(:,3,t+1).*spath(:,2,t+1)));
        invest = shrBore(t+1) - shrBore(t);
        highCost = norminv(shrBore(t+1)*P.inTruncProb+P.probBelow,P.investCostMean,P.investCostSD);
        if highCost<=lowCost
            investCost = 0;
        else
            investCost = P.investCostMean*invest+P.investCostSD^2/P.inTruncProb*(normpdf(lowCost,P.investCostMean,P.investCostSD)-normpdf(highCost,P.investCostMean,P.investCostSD));
        end

        gwUse(t,1) = (1-shrBore(t))*gwDug + shrBore(t).*gwBore;
        levelPath(t+1,1) = updateLevels(levelPath(t),gwUse(t),P);

        val(t,:) = [nbDug nbBore (1-shrBore(t))*nbDug + shrBore(t)*nbBore - investCost invest*P.convertTax];
        xPath(t,:) = [shrBore(t+1)-shrBore(t) gwDug gwBore];

        valChange = abs(val(t,3)*P.discount^(t-1));
        t = t+1;
        if t>=nyrs+1
            display('Stopping reSolve because we''ve hit max yrs instead of valTol')
            break
        end

    end

    compareLength = min(numel(assumedPath),numel(levelPath));
    plot(assumedPath)
    hold on;
    plot(levelPath)
    pause(1)
    pathChange = max(abs(assumedPath(1:compareLength)-levelPath(1:compareLength)));
    fprintf ('%4i %10.1e\n',iter,pathChange)
    
    cycleCount = cycleCount+1;
    if cycleCount>20
        moveFrac = moveFrac/2;
        cycleCount = 1;
    end
    assumedPath = (1-moveFrac)*assumedPath(1:compareLength) +moveFrac*levelPath(1:compareLength);
    

    %find maxCost of farms that have adopted by end of period
end

reOutput.statePath = [shrBore' levelPath];
reOutput.controlPath = xPath;
reOutput.valPath = val;
reOutput.reVal = (P.discount.^(0:length(val)-1))*val(:,3);
