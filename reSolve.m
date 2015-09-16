function reOutput = reSolve(P,modelOpts,paramGuess)
    if ~exist('paramGuess')
        projectedLevelsIn = P.h0*ones(modelOpts.trendPts,1);
        YrsLeft = 0;
        levelParams(modelOpts.trendPts) = P.h0;
    else
        YrsLeft = paramGuess.yrsLeft;
        levelParams = paramGuess.levelParams;
    end
    change = 100;
	model.func = 'cpFunc';
	model.discount = P.discount;
    while change>modelOpts.ttol
	    P.levelParams = levelParams;
        model.params = {P};
        smin = [0 0]; %both well capital and time remaining have natural minimums at 0
		smax = [P.maxWellCap max(YrsLeft,modelOpts.minT)]; %we may have trouble with the maxWellCap since there is no natural max
		%create the matrix of nodes
        n = [modelOpts.capNodes min(modelOpts.yrNodes,smax(2))];
		fspace = fundefn('cheb',n,smin,smax);
        snodes = funnode(fspace);
        s = gridmake(snodes);
        v = ones(size(s,1),1); x = [v v];
        [~,s,~,x] = dpsolve(model,fspace,s,v,x);
		[ssim,xsim] = dpsimul(model,[P.maxWellCap YrsLeft],YrsLeft,s,x);
        levelPath = P.h0;
        if YrsLeft
            for t=1:YrsLeft
                nextState = optFunc('g',[ssim(:,1,t) levelPath(t)],xsim(1,:,t),[],P);
                nextLevel(t) = nextState(2);
                levelPath(t+1) = nextLevel(t);
            end
        else
            nextLevel(1) = levelPath(1);
            levelPath(2) = levelPath(1);
        end
        oldLevels = polyval(levelParams,1:YrsLeft+1);
        change = max(abs(oldLevels - nextLevel));
        if YrsLeft>5
            rowInds = [1 2 floor(YrsLeft/2) YrsLeft];
        else
            rowInds = 1:YrsLeft;
        end
        if any(rowInds)
            matchYrs = rowInds-1;
            matchLevels = levelPath(rowInds);
        else
            matchYrs = 0:modelOpts.trendPts-1;
            matchLevels = nextLevel*ones(1,modelOpts.trendPts); 
        end
        [YrsLeft,levelParams] = computeTrend(matchLevels,matchYrs);        
    end
	
    %check whether our final action seems like we're really at a steady
    %state
    lastChange = levelPath(end) - levelPath(end-1)
    
    reOutput.cvValue = optFunc('f',[ssim(:,1) levelPath],xsim,[],P);
    reOutput.pvValue = (model.discount.^(0:YrsLeft+1)).*reOutput.cvValue;
    reOutput.capitalPath = ssim(:,1);
    reOutput.levelPath = levelPath;
	
