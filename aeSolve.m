function aeOutput = aeSolve(P,modelOpts)
	model.func = 'cpFunc';
	model.discount = P.discount;
    recentLevels = [P.h0-P.levelTrend; P.h0]; recentYears = [-1; 0];
    if numel(recentLevels)>2
        [YrsLeft,levelParams]=computeTrend(recentLevels,recentYears);
    else
        YrsLeft = (P.bottom-P.h0)/P.levelTrend;
        levelParams = [-P.levelTrend P.bottom];
    end
    aeVal = 0;
    npvValue = 0;
    t=0;
    
    statePath(P.wellInd) = P.wellCap0;
    statePath(P.levelInd) = P.h0;
    statePath(P.yrsInd) = YrsLeft;
    paramPath = zeros(1,modelOpts.trendPts)
    paramPath(1,1:length(levelParams))=levelParams;
    
	while npvValue>modelOpts.vtol||t<modelOpts.minT
        %compute the path that farmers will take as given based on recent
        %history
        P.levelParams = levelParams;
        model.params = {P};
        smin = [0 0]; %both well capital and time remaining have natural minimums at 0
        smax = [P.maxWellCap max(YrsLeft,modelOpts.minT)]; %we may have trouble with the maxWellCap since there is no natural max
        %create the matrix of nodes
        n = [modelOpts.capNodes min(modelOpts.yrNodes,smax(2))];
        fspace = fundefn('spli',n,smin,smax);
        snodes = funnode(fspace);
        s = gridmake(snodes);
        
        %solve the dynamic programming model, taking future water levels as
        %deterministic
        v = ones(size(s,1),1); x(:,P.investInd) = P.maxInvest*v; x(:,P.gwInd) = v;
        [~,s,~,x] = dpsolve(model,fspace,s,v,x);
        pseudoState = statePath(t+1,[P.wellInd P.levelInd]);
        simulState(1,:) = statePath(t+1,[P.wellInd P.yrsInd]);
        [ssim,xsim] = dpsimul(model,simulState,1,s,x);
        npvValue = model.discount^t*optFunc('f',pseudoState,xsim(1,:,1),[],P);
        thisAction = xsim(1,:,1);
        nextState = optFunc('g',pseudoState,thisAction,[],P);
        nextLevel = nextState(1,2);

        aeVal = aeVal + npvValue;
        t = t+1;
        % update the information used to compute the trend. If we haven't
        % accumulated the max history yet, we'll add points. If we've
        % already computed the max history, we'll replace the oldest point
        if length(recentLevels)<modelOpts.trendPts
            recentLevels = [recentLevels; nextLevel];
            recentYears = [1-length(recentLevels):0]';
        else
            recentLevels = [recentLevels(2:end); nextLevel];
        end
        [YrsLeft,levelParams]=computeTrend(recentLevels,recentYears);
        controlPath(t,:) = thisAction;
        statePath(t+1,P.wellInd) = ssim(1,1,end);
        statePath(t+1,P.yrsInd) = YrsLeft;
        statePath(t+1,P.levelInd) = nextLevel;
        paramPath(t+1,1:length(levelParams))=levelParams;
	end
		
	
