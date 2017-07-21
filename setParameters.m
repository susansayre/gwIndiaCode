close all; dbstop if error
csCross = 0;

parameterMat = {
    %1=ParameterName    2=ParameterDescription      					3=ParameterUnits        %4=BaseValue		%Note source following % sign
    'discount'          'discount rate'                                         '%'             1/1.1;
    'chokePrice'        'max MB of water'                                       'kRs/kmc'       2; %1.3            %assume same choke price
    'mbAtMax'           'mb at maxWater'                                        'kRs/kmc/kmc'   1.13;
    'maxWater'          'max vol to fully irrigate all land'                    'kmc/ha'        40;
    'maxCapacity'       '% achievable by each tech when not depth constrained'  'kmc/ha'        [1 1 1];
    'maxLiftMaxCap'     'maximum lift to get full cap for tech'                 'm'             [10 30 60];
    'maxDepths'         'max depth of each technology'                          'm'             [25 50 75];
    'energyNeed'        'energy needed to pump vol up 1m'                       'kWh/m3/m'      1/367;
    'pumpEfficiency'    ''                                                      '%'             [.25 .25 .25];
    'kWhCost'           'energy cost per kWh per m'                             'kRs/kmc/m'     [3.3 3.3 3.3];
    'eCostShr'          'share of real energy cost paid by farmers'             '%'             0;
    'avgPerKwH'         'average per kWh cost paid by farmers'                  'Rs/kWh'        .5;
    'landHeight'        'Initial Surface Level'                                 'm'             0;
    'minLift'           ''                                                      'm'             2;
    'initialLift'       'Initial m of pumping'                                  'm'             2; %comparative statics on this
    'bottom'            'Aquifer bottom'                                        'm'             -100;
    'maxExploitRatio'   'ratio max consumptive use to inflow'                   '%'             3;
    'returnFlow'        'Share of applied water that percolates back'           '%'             .25;
    'annualDrop'        'current annual drop'                                   'm'             3;
    'currIrrigShare'    'current share irrigated'                               '%'             1;
    'initShares'        'initial shrs in tech 2 and 3'                          '%'             [1 0 0];
    'baseIC'            'minimum investment cost tech 1 to 2'                   'kRs/ha'        40;
    'maxIC'             'maximum investment cost tech 1 to 2'                   'kRs/ha'        80;
    'ic1cost'           'ratios'                                                '% of base'     [0 1 1.6];
    'ic2cost'           'cost of moving from tech 2'                            '% of base'     [0 .04 1];
    'ic3cost'           'cost of moving from tech 3'                            '% of base'     [0 0 .1];              
    'deltaX'            'deviation for numeric derivs'                          ''              1e-5;
    'FEcostSocShr'      'shr of the fixed energy costs that count in soc'       ''              0;
%    'logisticRate'      ''                                                      ''              10;
    'viableTechs'       'indices of technologies being used'                    ''              [1 2 3];
    'interactive'       ''                                                      ''              0;
    'checkAtInvest'     ''                                                      '%'             .01;
    'maxInvest'         ''                                                      ''              1;
    };

numParams = size(parameterMat,1);
%initialize parameter values
for ii=1:numParams
    thisParam = parameterMat{ii,1};
    eval(['P.' thisParam '= parameterMat{ii,4};']) 
    %sets the value of each parameter listed in the first column to its value in the 4th column;
end

%comparative statics can be performed on as many variables as desired
%different types of comparative statics are coded by different types and
%require different inputs in the last column
%types
%1 = straight levels -- input = values used
%2 = level range -- input = [min max numPts(>=2)
%3 = straight % -- input = percent of base value
%4 = % range -- input = [minShare maxShare numPts(>=2)


compStatParams = {
    %parameterName      %type       %inputs
%     'pumpEfficiency(3)'    1           [.25 .3];
%        'discount'           1          [1/1.1];
%       'discount'           1          [1/1.1 1/1.05];
%       'baseIC'             1          [30 40 50];
%      'maxIC'              1          [60 80 100];
%        'maxIC'             1           100;
%        'discount'          1           [1/1.05];
%        'midICyears'         1          [1.5 3 4.5];
%        'meanMult'           1          [.5 1 2];
%       'discount'           1          [1/1.1 1/1.05];
%        'chokePrice'         1          [2 1.65 1.3];
%        'mbAtMax'             3          [1 .75 .5];
%        'maxLiftMaxCap(1)'    1          [5 10 15 20];
%       'midICyears'          1          [.75 1 1.25];
%      'meanMult'           1           [.5 .75 1];
        'maxExploitRatio'    1           [2 3];
        'annualDrop'         1           2;
%       'annualDrop'         1             [1.5];
%      'ic3cost'            1           
%      'annualDrop'          1           [ 1 2 3];
    };

numCompStatParams = size(compStatParams,1);
for ii=1:numCompStatParams
    switch compStatParams{ii,2}
        case 1
            csArray{ii} = compStatParams{ii,3}';
            csInds{ii} =(1:numel(compStatParams{ii,3}))';
        case 2
            min = compStatParams{ii,3}(1);
            max = compStatParams{ii,3}(2);
            step = (max-min)/(compStatParams{ii,3}(3)-1);
            csArray{ii} = (min:step:max)';
            csInds{ii} = (1:numel(csArray{ii}))';
        case 3
            eval(['baseValue = P.' compStatParams{ii,1} ';'])
            csArray{ii} = (compStatParams{ii,3}*baseValue)';
            csInds{ii} = (1:numel(compStatParams{ii,3}))';
        case 4
            eval(['baseValue = P.' compStatParams{ii,1} ';'])
            min = compStatParams{ii,3}(1);
            max = compStatParams{ii,3}(2);
            step = (max-min)/(compStatParams{ii,3}(3)-1);
            csArray{ii} = (min:step:max)'*baseValue;
            csInds{ii} = (1:numel(csArray{ii}))';
        otherwise
            error(['I don''t know case ' compStatParams{ii,2}])
    end
end

if csCross
    csValues = gridmake(csArray);
    csIndMat = gridmake(csInds);
else
    csValues = [];
    csIndMat = [];
    for ii=1:numel(csArray)
        baseInd = zeros(1,numel(csArray));
        baseValue = zeros(1,numel(csArray));
        for jj=1:numel(csArray)
            if jj==ii, continue, end
            eval(['baseValue(:,jj) = P.' compStatParams{jj,1} ';'])
        end
        for jj=1:numel(csArray{ii});
            thisValue = baseValue;
            thisInd = baseInd;
            thisValue(:,ii) = csArray{ii}(jj);
            thisInd(:,ii) = jj;
            csValues = [csValues; thisValue]; 
            csIndMat = [csIndMat; thisInd];
        end
    end
end

[cases,csParams] = size(csValues);

%set model options that won't change across cs runs
modelOpts.heightNodes = 25; %number of approximation nodes for levels, save one for maxDepth Dug
modelOpts.capNodes = 15; %number of approximation nodes for well capital
modelOpts.icMnodes = 100;
%modelOpts.yrNodes = 10;
modelOpts.minT = 100; %minimum number of years forward to simulate
modelOpts.vtol = 1e-2; %value function convergence tolerance
modelOpts.ttol = .5; %level trend convergence tolerance
modelOpts.stol = 1e-3; %invest path trend convergence tolerance
modelOpts.algorithm = 'funcit';
modelOpts.maxit = 300;
modelOpts.nres = 2;
modelOpts.refineExisting = 0;

modelOptsNi = modelOpts;
% modelOptsNi.capNodes = 4; %since investment won't happen, we don't need very many nodes
% modelOptsNi.icNodes = 4;
%Loop through compStat cases
if ~exist('runID','var')
    runID = datestr(now,'yyyymmdd_HHMMSS');
    doRun = 'Y';
else
    doRun = input(['Do you want to continue the existing run stored in ' runID '? Y/N [N]'],'s')
    if isempty(doRun)
        doRun = 'N';
    end
end

if ~strcmp(doRun,'Y')
    disp('Aborting run so I don''t overwrite results')
    return
end

for ii=1:cases
        thisP = P;
        thisP.csString = [];
        for jj=1:numCompStatParams
            %set parameter values based on compStat Case
            eval(['thisP.' compStatParams{jj,1} '= csValues(ii,jj);'])
            thisP.csString = [thisP.csString compStatParams{jj,1} ' = ' num2str(csValues(ii,jj))];
            if jj<numCompStatParams
                thisP.csString = [thisP.csString ', '];
            end
        end
    paramCases{ii} = thisP;
end

if ~exist('detailedOutput','dir')
    mkdir('detailedOutput')
end
if ~exist(['detailedOutput/' runID],'dir')
    mkdir(['detailedOutput/' runID])
end 
save(['detailedOutput/' runID '/setUp'])

%runTypes = {'', 'niP0', 'niPEnd', 'niPAvg', 'niPMid'};
runTypes = {'','startPt'};
runTypes = {'','only2','only1'};
runTypes = {''};

for ii=1:cases
    
    for jj = 1:numel(runTypes)
        runType = runTypes{jj};
        thisID = ['case' num2str(ii) runTypes{jj}];
        if ~modelOpts.refineExisting && exist(['detailedOutput/' runID '/' thisID '.mat'],'file')
            display([ thisID ' already computed. Loading results'])
            load(['detailedOutput/' runID '/' thisID],'output')
            eval(['results' runType '{ii} = output;'])
            eval(['pgainRE' runType '(ii) = output.reCp.pgain(1);'])
            clear output;
        else
            theseOpts = modelOptsNi; %will be overwritten for general case
            switch runType
                case ''
                    thisP = paramCases{ii};
                    theseOpts = modelOpts;
                 case 'only2'
                    thisP = paramCases{ii};
                    thisP.viableTechs = [1 2];
                    thisP.initShares = results{ii}.reRealE.sharePath(end,:);
                case 'only1'
                    thisP = paramCases{ii};
                    thisP.viableTechs = [1 2];
                    thisP.maxInvest = 0;
                case 'startPt'
                    thisP = paramCases{ii};
                    thisP.initialLift = thisP.maxLiftMaxCap(1);
                otherwise
                    warning(['I don''t recognize run type ' runType{jj} ])
                    continue
            end
            disp(['Solving ' thisID])
            eval(['results' runType '{ii} = solveCase(thisP,theseOpts,{runID thisID});'])
%             eval(['pgainRE' runType '(ii) = results' runType '{ii}.reCp.pgain(1)'])   
            close all
        end
    end
end
% diary(['detailedOutput/' runID '/summaryLog.txt'])

% compStatParams
% csValues
% %pgainAE
% pgainRE
% optVal = extractResult(results,'opt.val');
% optVal2 = extractResult(results,'opt.optVal')
% %aeVal = extractResult(results,'aeOut.aeVal')
% reVal = extractResult(results,'reCp.reVal')

% try
% %    aeValPath = extractResult(results,'aeOut.valPath(1:50,:)')
%     reValPath = extractResult(results,'reCp.valPath(1:50,:)')
%     optValPath = extractResult(results,'opt.valPath(1:50,:)')
%     
% %    aeStatePath = extractResult(results,'aeOut.statePath(1:50,:)')
%     reStatePath = extractResult(results,'reCp.statePath(1:50,:)')
%     optStatePath = extractResult(results,'opt.statePath(1:50,:)')
%     
% %    aecontrolPath = extractResult(results,'aeOut.controlPath(1:50,:)')
%     recontrolPath = extractResult(results,'reCp.controlPath(1:50,:)')
%     optcontrolPath = extractResult(results,'opt.controlPath(1:50,:)')   
% catch
% end

% diary off

save(['detailedOutput/' runID '/fullResults'])
