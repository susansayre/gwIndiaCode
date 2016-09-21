close all; dbstop if error
initializeClock = tic;
% if isempty(gcp)
%    % parpool
% end
initializeTime = toc(initializeClock);
programTime = tic;


parameterMat = {
    %1=ParameterName    2=ParameterDescription      					3=ParameterUnits        %4=BaseValue		%Note source following % sign
    'discount'          'discount factor (npv = discount^t*cv)'                 '%'             .95;
    'idDugInt'          'max MB of water'                                       'kRs/kmc'      .3*1.0357;          %2.5 k(m^3) increases rev from 7.5 kRs to 10 kRs time 30% profit from gross rev
    'idDugSlope'        ''                                                      'kRs/kmc/kmc'   .3*.0429;           %fit a curve to lots of point ranging from 0 to 5 k m^3
    'dugMax'            'max vol with dug well that doesn''t go dry'            'kmc/ha'        2.5;
    'idBoreInt'          'max MB of water'                                       'kRs/kmc'      .3*0.9294;          %2.5 k(m^3) increases rev from 7.5 kRs to 10 kRs time 30% profit from gross rev
    'idBoreSlope'        ''                                                      'kRs/kmc/kmc'   .3*.0121;          %fit a curve to lots of points ranging from 0 to 15 k m^3
    'boreMax'            'max vol with dug well that doesn''t go dry'            'kmc/ha'        15;                %this max is an artificial way of setting the MB to zero above this level 
    'electricityBore'   'cost per vol per m'                                    'kRs/kmc/m'     .012;
    'electricityDug'    ''                                                      'kRs/kmc/m'       0;
    'boreOM'            ''                                                      'kRs/ha'         .5;
    'landHeight'        'Initial Surface Level'                                 'm'             15;
    'minLift'           ''                                                      'm'             2;
    'initialLift'       'Initial m of pumping'                                  'm'             2; %comparative statics on this
    'bottom'            'Aquifer bottom'                                        'm'             -100;
    'inflow'            'Natural inflow per year'                               'kmc/ha/yr'      8;
    'returnFlow'        'Share of applied water that percolates back'           '%'             .35;
    'AS'                'area*storativity (vol released from 1 m of aquifer)'   'kmc/m/ha'       2.6;
    'shrBore0'          'initial shr in bore wells'                             '%'             .05;
    'maxDepthDug'       'max depth for a dug well'                              'm'             8;              %paper about poverty
    'investCostBase'    'common investment cost'                                'kRs/ha'        22.5;            %from Vis
    'investCostSD'      'standarized SD of indiv investment cost'               '$'             20;             %guess
    'investCostMean'    'mean of indiv investment cost'                         '$'             50;
    'convertTax'        'tax paid for drilling a well/converting'               '$/parcel'      0;   
    'slopeBoreZero'     'slope of the penalty for approaching bottom at zero'   ''              0;
    'minInvestCost'     'minimum value in the investCost distribution'          '$/parcel'      -Inf;
    'maxInvestCost'     'maximum value in the investCost distribution'          '$/parcel'      Inf;
    'investCostPenalty' 'penalty paid for high investment in one period'        '$/shr/shr'     0;
    'liftFullD'         'lift below which dug wells never go dry'               'm'             5; 
    'icDecayRate'       'rate at which common inv cost decays'                  '%'             .01;
    'boreLimitDecline'  'rate at which boreMax declines below liftFullD'        'k m^3/m'       0;
    'investLimit'       'upper bound on investment possible in one year'        '%'             1;
    'maxShr'            'maximum share in modern allowed'                       '%'             .99; %used to prevent problems with infinity
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
%    'electrictyDug'        1        [0 .05 .1 .15];
%    'investCostMean'       1        [5 10 15];
%     'inflow'               3        [2 3 1];
%     'inflow'               3        [1 2 3 4];
     'AS'                   3           [1 2];
%     'icDecayRate'          1        [.1 .07];
%     'investCostSD'         1        [10 8];
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

csValues = gridmake(csArray);
csIndMat = gridmake(csInds);

[cases,csParams] = size(csValues);

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
%Loop through compStat cases
if ~exist('runID','var')
    runID = datestr(now,'yyyymmdd_HHMMSS');
end

for ii=1:cases
        thisP = P;
        for jj=1:numCompStatParams
        %set parameter values based on compStat Case
        eval(['thisP.' compStatParams{jj,1} '= csValues(ii,jj);'])
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

runTypes = {'', 'niP0', 'niPEnd', 'niPAvg', 'niPMid'};

for ii=1:cases
    
    for jj = 1:numel(runTypes)
        runType = runTypes{jj};
        thisID = ['case' num2str(ii) runTypes{jj}];
        if exist(['detailedOutput/' runID '/' thisID '.mat'],'file')
            display([ thisID ' already computed. Loading results'])
            load(['detailedOutput/' runID '/' thisID],'output')
            eval(['results' runType '{ii} = output;'])
            eval(['pgainRE' runType '(ii) = output.reOut.pgain;'])
            clear output;
        else
            theseOpts = modelOptsNi; %will be overwritten for general case
            switch runType
                case ''
                    thisP = paramCases{ii};
                    theseOpts = modelOpts;
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
                    thisP.boreMax = wgtAvgShrBore*thisP.boreMax + (1-wgtAvgShrBore)*thisP.dugMax;
                case 'niPMid'
                    thisP = results{ii}.P;
                    shrBore = results{ii}.reOut.statePath(:,1);
                    wgtAvgShrBore = (thisP.discount.^(1:length(shrBore)))*shrBore/sum(thisP.discount.^(1:length(shrBore)));
                    thisP.investLimit = 0;
                    thisP.shrBore0 = wgtAvgShrBore;
                otherwise
                    warning(['I don''t recognize run type ' runType{jj} ])
                    continue
            end
            disp(['Solving ' thisID])
            eval(['results' runType '{ii} = solveCase(thisP,theseOpts,{runID thisID});'])
            eval(['pgainRE' runType '(ii) = results' runType '{ii}.reOut.pgain'])   
            close all
        end
    end
end
diary(['detailedOutput/' runID '/summaryLog.txt'])

compStatParams
csValues
%pgainAE
pgainRE
optVal = extractResult(results,'opt.val');
optVal2 = extractResult(results,'opt.optVal')
%aeVal = extractResult(results,'aeOut.aeVal')
reVal = extractResult(results,'reOut.reVal')

% try
% %    aeValPath = extractResult(results,'aeOut.valPath(1:50,:)')
%     reValPath = extractResult(results,'reOut.valPath(1:50,:)')
%     optValPath = extractResult(results,'opt.valPath(1:50,:)')
%     
% %    aeStatePath = extractResult(results,'aeOut.statePath(1:50,:)')
%     reStatePath = extractResult(results,'reOut.statePath(1:50,:)')
%     optStatePath = extractResult(results,'opt.statePath(1:50,:)')
%     
% %    aecontrolPath = extractResult(results,'aeOut.controlPath(1:50,:)')
%     recontrolPath = extractResult(results,'reOut.controlPath(1:50,:)')
%     optcontrolPath = extractResult(results,'opt.controlPath(1:50,:)')   
% catch
% end

diary off

save(['detailedOutput/' runID '/fullResults'])
timeToComplete = toc(programTime)
initializeTime