function reOutput = reSolve(P,modelOpts,initialPath)

 model.func = 'cpFunc';
 model.discount = P.discount;
 model.actions = (1:P.numTech)'; %convert to tech 1, convert to tech 2, convert to tech 3
 model.discretestates = [P.ind.tech P.ind.yrInd];
 model.e = 0; %shock values
 model.w = 1; %shock probabilities

P.numActions = numel(model.actions);
npvTol = modelOpts.vtol;
trendTol = modelOpts.ttol;
shrTol = modelOpts.stol;

levelPath = initialPath(1:end-1);

levelPathChange = trendTol*10;
iter = 0;

cycleCount = 1;
moveFrac = .2;
optset('dpsolve','showiters',0);
reOutput.converged = 0;
for iter=1:1000
    
    %solve common property dynamic optimization with level given
    P.levelPath = levelPath;
    P.lastYear = numel(levelPath);
    model.params = {P};
    
    fspace = fundefn('spli',modelOpts.icMnodes,0,P.maxICM,3,(1:P.numTech)',[1:length(levelPath)]');
    scoord = funnode(fspace);
    snodes = gridmake(scoord);

    s0(:,P.ind.icM) = scoord{P.ind.icM};

    tech0 = find(P.initShares);
    if numel(tech0)>1
        ub = 1;
        for ii=1:numel(tech0)
            lb = ub-P.initShares(tech0(ii));
            pickThese = intersect(find((1:modelOpts.icMnodes)./modelOpts.icMnodes>=lb),find((1:modelOpts.icMnodes)./modelOpts.icMnodes<=ub));
            s0(pickThese,P.ind.tech) = tech0(ii);
            ub = lb;
        end
    else
        s0(:,P.ind.tech) = tech0;
    end
    s0(:,P.ind.yrInd) = 1;
    
    if P.maxInvest
        %investment is allowed so we have to solve dynamic optimization and simulate forward
        vGuess = waterLimit(P.landHeight - levelPath(snodes(:,P.ind.yrInd)),P,snodes(:,P.ind.tech));
        [c,s,v,x] = dpsolve(model,fspace,snodes,vGuess);
        [ssim,xsim] = dpsimul(model,s0,P.lastYear-1,s,x);

        for ii=1:P.numTech
            sharePath(:,ii) = sum(squeeze(ssim(:,P.ind.tech,:)==ii))/length(s0);
        end
    else
        %no investment is allowed so we want to generate "simulation" results that match this outcome
        x = snodes(:,P.ind.tech);
        ssim(:,P.ind.icM,:) = repmat(s0(:,P.ind.icM),1,P.lastYear);
        ssim(:,P.ind.tech,:) = 1;
        ssim(:,P.ind.yrInd,:) = repmat(1:P.lastYear,size(ssim,1),1);
        xsim = repmat(s0(:,P.ind.tech),[1 1 P.lastYear]);   
        
        sharePath = repmat(P.initShares,P.lastYear,1);
    end
    levelPath = levelPath(1);
    for tt=1:P.lastYear
        lift = P.landHeight - levelPath(tt);

        wCosts = P.eCostShr*lift*P.eCosts; %ns x ntech

        maxWaters = waterLimit(lift,P);

        wCosts = P.eCostShr*lift*P.eCosts; %ns x ntech
    
        water(tt,:) = min(maxWaters,max((P.idInts - wCosts)./P.idSlopes,0));
            
        gw(tt,:) = sharePath(tt,:)*water(tt,:)';
        
        if tt<P.lastYear
            levelPath(tt+1,:) = updateLevels(levelPath(tt),gw(tt),P);
        end
    end
     
    levelPathChange = max(abs(P.levelPath-levelPath));
    
    fprintf ('%4i %10.1e\n',iter,levelPathChange)
%     figure()
%     plot(levelPath)
%     hold on;
%     plot(P.levelPath,'--')
%     pause(.001)
    
    levelPath = (1-moveFrac)*P.levelPath + moveFrac*levelPath; 
    
    if levelPathChange<trendTol
        reOutput.converged = 1;
        break;
    end
    iter = iter+1;
    %keyboard
    close
    
    cycleCount = cycleCount +1;
    if cycleCount>50
        moveFrac = moveFrac/2;
        cycleCount = 1;
    end
    
%    if iter>1000; keyboard; end
end

%trace benefits and costs over time for each farm
yrs = size(ssim,3); farms = size(ssim,1);
ssimLong = reshape(permute(ssim,[1 3 2]),yrs*farms,3);
xsimLong = reshape(permute(xsim,[1 3 2]),yrs*farms,1);
valPaths = reshape(cpFunc('f',ssimLong,xsimLong,[],P),farms,yrs);
npvs = sum(valPaths.*repmat(P.discount.^(0:1:yrs-1),farms,1),2);
valPriv = mean(valPaths)';
npvPriv = mean(npvs);
%generate states and actions that can be used to call optFunc to generate
%social costs and benefits

sharesRaw = share2RawShr(sharePath);

sPath(:,P.ind.level) = levelPath;
sPath(:,P.shareInds) = sharesRaw;
moveShr = zeros(P.numTech,P.numTech-1,yrs);
    
investShr = moveShr;
for ii=1:P.numTech
    posShrs = find(sharePath(:,ii));
    for jj=1:P.numTech-1
        moveShr(ii,jj,posShrs) = sum(squeeze(ssim(:,P.ind.tech,posShrs)==ii).*(xsim(:,posShrs)==jj))'./(farms*sharePath(posShrs,ii));
        if jj==1
            denominator = ones(numel(posShrs),1);
            doThese = posShrs;
        else
            denominator = squeeze(prod(1 - investShr(ii,1:jj-1,posShrs),2));
            doThese = posShrs(find(denominator));
        end
        investShr(ii,jj,doThese) = squeeze(moveShr(ii,jj,doThese))./denominator(find(denominator));
    end
end 
investPath = reshape(permute(investShr,[3 2 1]),yrs,P.numTech*(P.numTech-1));

xPath(:,P.ind.water) = water;
xPath(:,P.ind.invest) = investPath;

reOutput.statePath = sPath;
reOutput.controlPath = xPath;
reOutput.sharePath = sharePath;
reOutput.levelPath = levelPath;
reOutput.waterPath = water;

[nbS,nbDetailS,farmDetailS] = optValue(sPath,xPath,[],P,'society');
[nbI,nbDetailI,farmDetailI] = optValue(sPath,xPath,[],P,'indiv');
reOutput.nbPath = [nbS nbI];
reOutput.nbDetailPath = [nbDetailS nbDetailI];
reOutput.farmDetailS = farmDetailS;
reOutput.farmDetailI = farmDetailI;
reOutput.npv = P.discount.^(0:length(reOutput.nbPath)-1)*reOutput.nbPath;

%check for consistency
% [a,b] = sort(xsim,1,'descend'); %identifies which farms are selecting which technologies. Our assumptions imply the cheapest farms should be the ones in tech 3, mid in tech 2, and highest in tech1
% orderSwitches = find(b-repmat((1:size(b,1))',1,size(b,2)));
% if any(orderSwitches)
%     warning('The technology and cost order is not remaining constant')
% end

%plot policy functions
myStates.icM = snodes(:,P.ind.icM);
myStates.tech = snodes(:,P.ind.tech);
myStates.lift = P.h0 - levelPath(snodes(:,P.ind.yrInd));
myStates.choice = reshape(x,size(myStates.icM));
choicePlot = gramm('x',myStates.lift,'y',myStates.icM,'color',myStates.choice);
choicePlot.geom_point();
choicePlot.facet_grid(myStates.tech,[]);
choicePlot.set_names('x','water level','y','cost','color','choice','row','current tech');
choicePlot.set_color_options('map','brewer1');
choicePlot.set_point_options('base_size',10);

reOutput.choicePlot = choicePlot;