parameterMat = {
    %1=ParameterName    2=ParameterDescription      					3=ParameterUnits    %4=BaseValue		%Note source following % sign
    'wellCapDecayRate'  '% of well capital that must be replaced each year'     '%'             .01 ; %assumption
    'maxWellCap'        'upper bound on well capital'                           '%'             1;
    'discount'          'discount factor (npv = discount^t*cv)'                 '%'             .95;
    'gwdIntercept'      'Intercept of gw demand curve'                          '$/acre-inch'	10;
    'gwdSlope'          'Absolute value of slope of gw demand'                  '$/ai^2'    	1;
    'electricity'       'cost per hacm per m'                                   '$/m/hacm'      1;
    'maxDepthNoCap'     'maximum depth that a no capital well can reach'        'm'             8;
    'varInvestCost'     ''                                                      '$'             1;
    'fixedInvestCost'   ''                                                      '$'             0;
    'landHeight'        'Initial Surface Level'                                 'm'             500;
    'initialLift'       'Initial m of pumping'                                  'm'             5;
    'bottom'            'Aquifer bottom'                                        'm'             0;
    'inflow'            'Natural inflow per year'                               'ac-in/yr'      15.86*.08;
    'returnFlow'        'Share of applied water that percolates back'           '%'             .02;
    'AS'                'area*storativity (ha.cm released from 1 m of aquifer)' 'hacm/m'        10;
    'maxInvest'         'max investment possible in one year'                   '$'             .1;
    'levelTrend'        'linear level change rate/yr for use in AE model'       'm'             -1;
    'wellCap0'          'initial well capital level'                            '$'             0;
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
    'fixedInvestCost'   1           [0 1];
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
modelOpts.heightNodes = 5; %number of approximation nodes for levels
modelOpts.capNodes = 5; %number of approximation nodes for well capital
modelOpts.yrNodes = 5;
modelOpts.minT = 25; %minimum number of years forward to simulate
modelOpts.trendPts = 4; %maximum number of point used to compute trend expectations
modelOpts.vtol = 1; %value function convergence tolerance
modelOpts.ttol = 1; %level trend convergence tolerance
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