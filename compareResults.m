
compStatCases = { %date
    '20161024_094141'   1   'baseline';
    '20161102_135035'   1   'lower min cost';
    '20161102_135035'   3   'higher min cost';
    '20161102_135035'   4   'lower max cost';
    '20161102_135035'   6   'higher max cost';
%     '20161103_131554'   2   '5% discount rate';
% %    '20161103_145806'   1   '5% discount rate and higher max cost';
%     '20161026_145120'   3   'lower max marginal benefit';
%     '20161026_145120'   4   'lower marginal benefit at max';
%     '20161028_130227'   2   'lower annual drop';
%     '20161106_203516'   1   'higher recharge';
% %    '20161026_114817'   1   'lower maintenance costs';
};
    %not used
%    '20161022_212008'   'midICyears and meanMult'

numCompCases = size(compStatCases,1);

resultTypes = {
    'endLevels'     [numCompCases 3]      [-75 0]   'Water level at year 50'                                    'm';
    'shareFarming'  [numCompCases 3]      [0 1]     'Share irrigating at year 50'                               '%';
    'npvs'          [numCompCases 3]      [0 600]   'Net present value of irrigation benefits'                  'Rs./ha';
    'npvRevenue'    [numCompCases 3]      [0 1200]  'Present value of gross benefits'                           'Rs./ha';
    'npvEnergy'     [numCompCases 3]      [0 1000]  'Present value of pumping cost expenditures'                'Rs./ha';
    'npvInvest'     [numCompCases 3]      [0 200]   'Present value of investment and maintenance expenditures'  'Rs./ha';
    'pgain'         [numCompCases 2]      [0 1]    'Percentage gain from optimal management'                    '%';
    'endShares'     [numCompCases 3 3]    [0 1]     'Share using each technology at year 50'                    '%';
    'yrsToSS'       [numCompCases 3]      [0 100]   'Years to steady state depth'                               'yrs';
};


endYear = 60;

manageCases = {'reCp','reRealE','opt'};

for ii=1:numCompCases
    load(fullfile('detailedOutput',compStatCases{ii,1},['case' num2str(compStatCases{ii,2})]),'output','P')
    npvYears = 100;
    discountFactorRow = P.discount.^(0:1:npvYears-1);
    for jj=1:numel(manageCases)
        endLevels(ii,jj) = eval(['output.' manageCases{jj} '.levelPath(endYear,1)']);
        npvs(ii,jj) = eval(['output.' manageCases{jj} '.npv(1)']);
        if jj<numel(manageCases)
            pgain(ii,jj) = eval(['output.' manageCases{jj} '.pgain(1)']);
        else
            pgain(ii,jj) = 0;
        end
        endShares(ii,jj,:) = eval(['output.' manageCases{jj} '.sharePath(endYear,:)']);
        npvInvest(ii,jj) = eval(['discountFactorRow*output.' manageCases{jj} '.nbDetailPath(1:npvYears,4)']);
        npvRevenue(ii,jj) = eval(['discountFactorRow*output.' manageCases{jj} '.nbDetailPath(1:npvYears,1)']);
        npvEnergy(ii,jj) = eval(['discountFactorRow*output.' manageCases{jj} '.nbDetailPath(1:npvYears,2)']);
        
        if -endLevels(ii,jj)>P.maxDepths(1)
            if jj==3; keyboard; end;
            endShares(ii,jj,1) = 0;
        end
        
        levelChange(ii,jj,:) = eval(['output.' manageCases{jj} '.levelPath(2:endYear) - output.' manageCases{jj} '.levelPath(1:endYear-1)']);
    end
end

levelChange = abs(levelChange)>.01;
yrsToSS = sum(levelChange,3);
shareFarming = sum(endShares,3);

myColors = [ 141 160 203; 252 141 98; 102 194 165]/255;
combinedFig = figure();
set(gcf,'PaperPosition',[1 1 6.5 9])
for ii=1:6
    subplot(4,2,ii)
    thisPlot = eval(['bar(' resultTypes{ii,1} ')']);
    for jj=1:numel(thisPlot)
        set(thisPlot(jj),'FaceColor',myColors(jj,:));
    end
    ylabel(resultTypes{ii,5});
    hold on;
    eval(['plot(repmat(xlim'',1,3),repmat(' resultTypes{ii,1} '(1,:),2,1),''k:'')'])
    set(gca,'YLim',resultTypes{ii,3});
    yLims = resultTypes{ii,3};
    set(gca,'YTick',[yLims(1) 0.5*(yLims(1)+yLims(2)) yLims(2)])
    title(resultTypes{ii,4})
    if ii<6
        set(gca,'XTickLabel',{})
    else
        set(gca,'XTickLabel',compStatCases(:,3))
        set(gca,'XTickLabelRotation',-90)
        legHand= legend('CP Flat Tariffs','CP Real Energy Costs','Optimal Managment');
        defaultPosition = get(legHand,'Position');
         defaultPosition(1:2) = [.15 .1];
         set(legHand,'Position',defaultPosition)
    end
end
saveas(gcf,'compStatBars','epsc')

% figure()
% plotBarStackGroups(endShares,compStatCases(:,3))

    
    
    