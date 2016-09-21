%fix reSolve problem
load('detailedOutput/20160810_144021/setUp')
runTypes = {'', 'niP0', 'niPEnd', 'niPAvg'};

%set model options that won't change across cs runs
modelOpts.heightNodes = 15; %number of approximation nodes for levels, save one for maxDepth Dug
modelOpts.capNodes = 15; %number of approximation nodes for well capital
modelOpts.icNodes = 10;
%modelOpts.yrNodes = 10;
modelOpts.minT = 100; %minimum number of years forward to simulate
modelOpts.vtol = 1e-2; %value function convergence tolerance
modelOpts.ttol = 1e-2; %level trend convergence tolerance
modelOpts.stol = 1e-3; %invest path trend convergence tolerance
modelOpts.algorithm = 'funcit';
modelOpts.maxit = 300;

modelOptsNi = modelOpts;
modelOptsNi.capNodes = 4; %since investment won't happen, we don't need very many nodes
modelOptsNi.icNodes = 4;

for ii=1:cases
    
    for jj = 1:numel(runTypes)
        runType = runTypes{jj};
        thisID = ['case' num2str(ii) runTypes{jj}];
        theseOpts = modelOptsNi; 
        doFullRun = 1;%will be overwritten for general case
        switch runType
            case ''
                doFullRun = 0;
            case 'niP0'
                thisP = results{ii}.P;
                thisP.investLimit = 0;
            case 'niPEnd'
                thisP = results{ii}.P;
                thisP.investLimit = 0;
                thisP.shrBore0 = results{ii}.reOut.statePath(end,1);
            case 'niPAvg'
                thisP = results{ii}.P;
                thisP.shrBore0 = 1; %all farms set to a wgt average farm type
                shrBore = results{ii}.reOut.statePath(:,1);
                wgtAvgShrBore = (thisP.discount.^(1:length(shrBore)))*shrBore/sum(thisP.discount.^(1:length(shrBore)));
                thisP.boreLimitDecline = (1-wgtAvgShrBore)*thisP.dugMax/(thisP.maxDepthDug - thisP.liftFullD);
                thisP.idBoreInt = wgtAvgShrBore*thisP.idBoreInt + (1-wgtAvgShrBore)*thisP.idDugInt;
                thisP.idBoreSlope = wgtAvgShrBore*thisP.idBoreSlope + (1-wgtAvgShrBore)*thisP.idDugSlope;
                thisP.electricityBore = wgtAvgShrBore*thisP.electricityBore + (1-wgtAvgShrBore)*thisP.electricityDug;
                thisP.investLimit = 0;
                thisP.maxShr = 1;
            otherwise
                warning(['I don''t recognize run type ' runType{jj} ])
                continue
        end
        disp(['Solving ' thisID])
        if doFullRun
            eval(['results' runType '{ii} = solveCase(thisP,theseOpts,{runID thisID});'])
        else
            results{ii} = redoBaseCase(runID,thisID);
        end
        eval(['pgainRE' runType '(ii) = results' runType '{ii}.reOut.pgain'])   
    end
end