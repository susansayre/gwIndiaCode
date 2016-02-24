function reOutput = reSolve(P,modelOpts,pathGuess,controlGuess)
    if ~exist('pathGuess')
        P.levelPath(1,:) = P.h0*ones(modelOpts.minT);
    else
        P.levelPath(1,:) = pathGuess;
    end

    P.explicitPath = 1;
    Yrs = length(P.levelPath);
    change = 100;
	model.func = 'cpFunc';
	model.discount = P.discount;
    iter = 0;
   % optset('dpsolve','algorithm','funcit');
    while change>modelOpts.ttol
        iter = iter + 1;
        smin = [0 0]; %both well capital and time remaining have natural minimums at 0
		smax = [1 Yrs]; 
		%create the matrix of nodes
        n = [modelOpts.capNodes modelOpts.yrNodes];
        n = [10 10];
		fspace = fundefn('spli',n,smin,smax);
        %stack an extra break point at 0 years left.
%         keyboard
        fspace.parms{2}{1} = sort([fspace.parms{2}{1}; 0]); %won't work unless I make a different guess for the value function at the second node.
        snodes = funnode(fspace);
        s = gridmake(snodes);
        model.params = {P};
        
        yrInds = floor(s(:,2));
        
        if exist('controlGuess')
            yrInds = max(1,length(controlGuess)-yrInds); % set the indicator to either yrs left or the max yrs left available in the control guess
            xGuess = controlGuess(yrInds,:);
        end
        if ~exist('v'); v = ones(size(s,1),1); v(1:length(snodes{1}))=0; end;
        if ~exist('xGuess')
            if ~exist('controlGuess')
                xGuess(:,P.investInd) = P.maxInvest*v; xGuess(:,P.gwDugInd) = v; xGuess(:,P.gwBoreInd) = v;
            else
                xGuess = zeros(size(s,1),3);
                evenYears = floor(s(:,2));
                xGuess(find(evenYears>0)) = controlGuess(length(controlGuess)+2-evenYears(find(evenYears>0)));
            end
        end
                
        [~,sOut,vOut,x] = dpsolve(model,fspace,s,v,xGuess);
        vOut = v;
		[ssim,xsim] = dpsimul(model,[P.shrBore0 Yrs],Yrs+5,sOut,x);
        levelPath = P.h0;
        for t=1:Yrs+5
            nextState = optFunc('g',[ssim(:,1,t) levelPath(t)],xsim(1,:,t),[],P);
            nextLevel(t) = nextState(2);
            levelPath(t+1) = nextLevel(t);
            if nextLevel(t)==0; break; end
        end
        oldYrs=Yrs;
        Yrs = length(levelPath);
        xGuess = xsim;
        compareInds = min(oldYrs,Yrs);
        change = max(abs(P.levelPath(1:compareInds)-levelPath(1:compareInds)));
        P.levelPath = levelPath;
    end
	
    %%I'm trying to improve the guesses, but it's not working.
    keyboard
    %check whether our final action seems like we're really at a steady
    %state
    lastChange = levelPath(end) - levelPath(end-1)
    
    reOutput.cvValue = optFunc('f',[squeeze(ssim(:,1,:))' levelPath'],xsim,[],P);
    reOutput.pvValue = (model.discount.^(0:YrsLeft+1)).*reOutput.cvValue;
    reOutput.capitalPath = ssim(:,1);
    reOutput.levelPath = levelPath;
	
