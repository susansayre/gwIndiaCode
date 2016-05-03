clear all; close all; dbstop if error
parameterMat = {
    %1=ParameterName    2=ParameterDescription      					3=ParameterUnits        %4=BaseValue		%Note source following % sign
    'discount'          'discount factor (npv = discount^t*cv)'                 '%'             .95;
    'dDugInt'           'Intercept of gw demand curve for dug wells'            '$/vol'         2;
    'dDugSlope'         'Absolute value of slope of gw demand for dug wells'    '$/vol^2'    	2;
    'dBoreInt'          'Intercept of gw demand curve for bore wells'            '$/vol'        1.8;
    'dBoreSlope'        'Absolute value of slope of gw demand for bore wells'    '$/vol^2'    	1;
    'electricity'       'cost per vol per m'                                    '$/m/vol'       .15;
    'landHeight'        'Initial Surface Level'                                 'm'             20;
    'initialLift'       'Initial m of pumping'                                  'm'             2;
    'bottom'            'Aquifer bottom'                                        'm'             0;
    'inflow'            'Natural inflow per year'                               'vol/yr'        .1;
    'returnFlow'        'Share of applied water that percolates back'           '%'             .02;
    'AS'                'area*storativity (vol released from 1 m of aquifer)'   'vol/m'         1;
    'maxInvest'         'max investment possible in one year'                   '%'             .1;
    'levelTrend'        'linear level change rate/yr for use in AE model'       'm'             -.1;
    'shrBore0'          'initial shr in bore wells'                             '%'             .01;
    'maxDepthDug'       'max depth for a dug well'                              'm'             8;              %paper about poverty
    'slopeMaxDepth'     'slope of cost function at max depth'                   '$/m'           2.5;         %completely made up
    'investCostMean'    'mean of investment cost distribution'                  '$/parcel'      3.0;            %guess
    'investCostSD'      'sd of investment cost'                                 '$'             1.5;             %guess
    'convertTax'        'tax paid for drilling a well/converting'               '$/parcel'      0;   
    'slopeBoreZero'     'slope of the penalty for approaching bottom at zero'   ''              0;
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
%    'AS'                1           [10   12.5    15];
   'dBoreInt'            1           [1.8];

%    'investCostMean'    1           [20 30 40];
%      'investCostMean'       1           [50 30 10];
%      'investCostSD'         1           [15 30 45];
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
modelOpts.heightNodes = 25; %number of approximation nodes for levels, save one for maxDepth Dug
modelOpts.capNodes = 25; %number of approximation nodes for well capital
%modelOpts.yrNodes = 10;
modelOpts.minT = 50; %minimum number of years forward to simulate
modelOpts.trendPts = 1; %maximum number of point used to compute trend expectations
modelOpts.vtol = 1e-1; %value function convergence tolerance
modelOpts.ttol = 1; %level trend convergence tolerance
modelOpts.algorithm = 'newton';
modelOpts.maxit = 2000;
modelOpts.icSteps = 9999;
%Loop through compStat cases
runID.timeStamp = datestr(now,'yyyymmdd_HHMMSS');
for ii=1:cases
    for jj=1:numCompStatParams
        %set parameter values based on compStat Case
        eval(['P.' compStatParams{jj,1} '= csValues(ii,jj);'])
    end
    %run problem for this parameter set
    runID.case = ['case' num2str(ii)];
    results{ii} = solveCase(P,modelOpts,runID);
    pgainAE(ii) = results{ii}.aeOut.pgain;
    pgainRE(ii) = results{ii}.reOut.pgain;
end
diary(['detailedOutput/' runID.timeStamp '/summaryLog.txt'])

compStatParams
csValues
pgainAE
pgainRE
optVal = extractResult(results,'opt.val')
optVal2 = extractResult(results,'opt.optVal')
aeVal = extractResult(results,'aeOut.aeVal')
reVal = extractResult(results,'reOut.reVal')

try
    aeValPath = extractResult(results,'aeOut.valPath(1:50,:)')
    reValPath = extractResult(results,'reOut.valPath(1:50,:)')
    optValPath = extractResult(results,'opt.valPath(1:50,:)')
    
    aeStatePath = extractResult(results,'aeOut.statePath(1:50,:)')
    reStatePath = extractResult(results,'reOut.statePath(1:50,:)')
    optStatePath = extractResult(results,'opt.statePath(1:50,:)')
    
    aecontrolPath = extractResult(results,'aeOut.controlPath(1:50,:)')
    recontrolPath = extractResult(results,'reOut.controlPath(1:50,:)')
    optcontrolPath = extractResult(results,'opt.controlPath(1:50,:)')   
catch
end

diary off

save(['detailedOutput/' runID.timeStamp '/fullResults'])