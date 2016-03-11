clear all; close all; dbstop if error
parameterMat = {
    %1=ParameterName    2=ParameterDescription      					3=ParameterUnits        %4=BaseValue		%Note source following % sign
    'discount'          'discount factor (npv = discount^t*cv)'                 '%'             .95;
    'dDugInt'           'Intercept of gw demand curve for dug wells'            '$/vol'         2;
    'dDugSlope'         'Absolute value of slope of gw demand for dug wells'    '$/vol^2'    	.2;
    'dBoreInt'          'Intercept of gw demand curve for dug wells'            '$/vol'         1.8;
    'dBoreSlope'        'Absolute value of slope of gw demand for dug wells'    '$/vol^2'    	.1;
    'electricity'       'cost per vol per m'                                    '$/m/vol'       .14;
    'landHeight'        'Initial Surface Level'                                 'm'             16;
    'initialLift'       'Initial m of pumping'                                  'm'             .1*16;
    'bottom'            'Aquifer bottom'                                        'm'             0;
    'inflow'            'Natural inflow per year'                               'vol/yr'        1;
    'returnFlow'        'Share of applied water that percolates back'           '%'             .02;
    'AS'                'area*storativity (vol released from 1 m of aquifer)'   'vol/m'         20;
    'maxInvest'         'max investment possible in one year'                   '%'             .1;
    'levelTrend'        'linear level change rate/yr for use in AE model'       'm'             -1.6;
    'shrBore0'          'initial shr in bore wells'                             '%'             .01;
    'maxDepthDug'       'max depth for a dug well'                              'm'             8;              %paper about poverty
    'slopeMaxDepth'     'slope of cost function at max depth'                   '$/m'           2.5;         %completely made up
    'investCostMean'    'mean of investment cost distribution'                  '$/parcel'      15;            %guess
    'investCostSD'      'sd of investment cost'                                 '$'             12;             %guess
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
    'initialLift'   3           [1 2];
    };

csArray = {};
csInds = {};
numCompStatParams = size(compStatParams,1);
for ii=1:numCompStatParams
    switch compStatParams{ii,2}
        case 1
            csArray = {csArray,compStatParams{ii,3}};
            csInds = {csInds,1:numel(compStatParams)};
        case 2
            min = compStatParams{ii,3}(1);
            max = compStatParams{ii,3}(2);
            step = (max-min)/(compStatParams{ii,3}(3)-1);
            csArray = {csArray,min:step:max};
            csInds = {csInds,1:compStatParams{ii,3}};
        case 3
            eval(['baseValue = P.' compStatParams{ii,1} ';'])
            csArray = {csArray,compStatParams{ii,3}*baseValue};
            csInds = {csInds,1:numel(compStatParams)};
        case 4
            eval(['baseValue = P.' compStatParams{ii,1} ';'])
            min = compStatParams{ii,3}(1);
            max = compStatParams{ii,3}(2);
            step = (max-min)/(compStatParams{ii,3}(3)-1);
            csArray = {csArray,(min:step:max)*baseValue};
            csInds = {csInds,1:compStatParams{ii,3}};
        otherwise
            error(['I don''t know case ' compStatParams{ii,2}])
    end
end

csValues = gridmake(csArray);
csIndMat = gridmake(csInds);

[cases,csParams] = size(csValues);

%set model options that won't change across cs runs
modelOpts.heightNodes = 10; %number of approximation nodes for levels, save one for maxDepth Dug
modelOpts.capNodes = 10; %number of approximation nodes for well capital
modelOpts.yrNodes = 10;
modelOpts.minT = 50; %minimum number of years forward to simulate
modelOpts.trendPts = 1; %maximum number of point used to compute trend expectations
modelOpts.vtol = 1; %value function convergence tolerance
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
end