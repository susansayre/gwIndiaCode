clear all; close all; dbstop if error
clear all; close all;
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
    'idDugSlope'        ''                                                      'kRs/kmc/kmc'   .3*.0429;
    'dugMax'            'max vol with dug well that doesn''t go dry'            'kmc/ha'        3;
    'idBoreInt'          'max MB of water'                                       'kRs/kmc'      .3*0.9294;          %2.5 k(m^3) increases rev from 7.5 kRs to 10 kRs time 30% profit from gross rev
    'idBoreSlope'        ''                                                      'kRs/kmc/kmc'   .3*.0121;
    'boreMax'            'max vol with dug well that doesn''t go dry'            'kmc/ha'        15;                %this max is an artificial way of setting the MB to zero above this level 
    'electricityBore'   'cost per vol per m'                                    'kRs/kmc/m'     .012;
    'electricityDug'    ''                                                      '$/m/vol'       0;
    'landHeight'        'Initial Surface Level'                                 'm'             15;
    'minLift'           ''                                                      'm'             .5;
    'initialLift'       'Initial m of pumping'                                  'm'             2;
    'bottom'            'Aquifer bottom'                                        'm'             -100;
    'inflow'            'Natural inflow per year'                               'kmc/ha/yr'        2;
    'returnFlow'        'Share of applied water that percolates back'           '%'             .35;
    'AS'                'area*storativity (vol released from 1 m of aquifer)'   'kmc/m/ha'       1.6;
    'maxInvest'         'max investment possible in one year'                   '%'             .1;
    'levelTrend'        'linear level change rate/yr for use in AE model'       'm'             0;
    'shrBore0'          'initial shr in bore wells'                             '%'             .05;
    'maxDepthDug'       'max depth for a dug well'                              'm'             8;              %paper about poverty
    'investCostMean'    'mean of investment cost distribution'                  'kRs/ha'        25;            %guess
    'investCostSD'      'standarized SD of investment cost'                     '$'             4.5;             %guess
    'convertTax'        'tax paid for drilling a well/converting'               '$/parcel'      0;   
    'slopeBoreZero'     'slope of the penalty for approaching bottom at zero'   ''              0;
    'minInvestCost'     'minimum value in the investCost distribution'          '$/parcel'      5;
    'maxInvestCost'     'maximum value in the investCost distribution'          '$/parcel'      40;
    'investCostPenalty' 'penalty paid for high investment in one period'        '$/shr/shr'     0;
    'depthFullD'        'depth above which dug wells never go dry'              'm'             5;
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
%    'AS'                3           [1];
%    'electrictyDug'        1        [0 .05 .1 .15];
%    'investCostMean'       1        [5 10 15];
    'inflow'               1        [2 4 6 8 10];
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
modelOpts.heightNodes = 30; %number of approximation nodes for levels, save one for maxDepth Dug
modelOpts.capNodes = 30; %number of approximation nodes for well capital
%modelOpts.yrNodes = 10;
modelOpts.minT = 100; %minimum number of years forward to simulate
modelOpts.trendPts = 1; %maximum number of point used to compute trend expectations
modelOpts.vtol = 1e-2; %value function convergence tolerance
modelOpts.ttol = 1e-2; %level trend convergence tolerance
modelOpts.algorithm = 'funcit';
modelOpts.maxit = 300;
modelOpts.maxitncp = 300;
%Loop through compStat cases
runID = datestr(now,'yyyymmdd_HHMMSS');
for ii=1:cases
    for jj=1:numCompStatParams
        %set parameter values based on compStat Case
        thisP = P;
        eval(['thisP.' compStatParams{jj,1} '= csValues(ii,jj);'])
    end
    paramCases{ii} = thisP;
end

for ii=1:cases
    %run problem for this parameter set
    thisID = ['case' num2str(ii)];
    results{ii} = solveCase(paramCases{ii},modelOpts,{runID thisID});
    %pgainAE(ii) = results{ii}.aeOut.pgain;
    pgainRE(ii) = results{ii}.reOut.pgain;
end
diary(['detailedOutput/' runID '/summaryLog.txt'])

compStatParams
csValues
%pgainAE
pgainRE
optVal = extractResult(results,'opt.val')
optVal2 = extractResult(results,'opt.optVal')
%aeVal = extractResult(results,'aeOut.aeVal')
reVal = extractResult(results,'reOut.reVal')

try
%    aeValPath = extractResult(results,'aeOut.valPath(1:50,:)')
    reValPath = extractResult(results,'reOut.valPath(1:50,:)')
    optValPath = extractResult(results,'opt.valPath(1:50,:)')
    
%    aeStatePath = extractResult(results,'aeOut.statePath(1:50,:)')
    reStatePath = extractResult(results,'reOut.statePath(1:50,:)')
    optStatePath = extractResult(results,'opt.statePath(1:50,:)')
    
%    aecontrolPath = extractResult(results,'aeOut.controlPath(1:50,:)')
    recontrolPath = extractResult(results,'reOut.controlPath(1:50,:)')
    optcontrolPath = extractResult(results,'opt.controlPath(1:50,:)')   
catch
end

diary off

save(['detailedOutput/' runID '/fullResults'])
timeToComplete = toc(programTime)
initializeTime