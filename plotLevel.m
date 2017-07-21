levelFig = figure();
set(levelFig,'PaperPosition',[1 1 3.25 5])
plotYears = 50;
levelData = squeeze(levels(1:plotYears,:,:));
allLines = levelData;
allLines(:,4) = -P.maxDepths(1);
allLines(:,5) = -P.maxDepths(2);
allLines(:,6) = -max(P.idInts)*max(P.pumpEfficiency)/(P.energyNeed*min(P.kWhCost));

%create water level graph
myColors = [27 158 119; 217 95 2; 117 112 179]/255;
myStyles = {'-','--','-.'};
%cheat by plotting the levels twice so that we get two y-axes
[levelPlot, h1, h2] = plotyy(1:plotYears,allLines,1:plotYears,levelData);
for ii=1:numel(h1)
    if ii<=numel(h2)
        set([h1(ii);h2(ii)],'Color',myColors(ii,:))
        set([h1(ii); h2(ii)],'LineWidth',1.5)
        set([h1(ii); h2(ii)],'LineStyle',myStyles{ii});
    else
        set(h1(ii),'Color',[0 0 0])
        set(h1(ii),'LineStyle',':')
        set(h1(ii),'LineWidth',.75)
    end
    set(levelPlot,'YLim',[-80 0])
end

myLegend = legend(flip(get(levelPlot(2),'Children')),'Optimal Management','CP Real Cost','CP Flat Tariff');
set(myLegend,'Location','SouthOutside')

%extract tick mark levels
tickPoints = get(levelPlot(1),'YTick');
for ii=1:numel(tickPoints)
    costLabels{ii} = num2str(P.eCosts(1)*(P.landHeight-tickPoints(ii)),'%1.1f');
end
set(levelPlot(2),'YTickLabel',costLabels)
set(get(levelPlot(2),'YLabel'),'String','Pumping Cost (Rs./m^{3})')
set(get(levelPlot(1),'YLabel'),'String','Water level (m below surface)')
xlabel('years')

text(25,allLines(1,4),'max depth tech 1','BackgroundColor',[1 1 1])
text(25,allLines(1,5),'max depth tech 2','BackgroundColor',[1 1 1])
text(25,allLines(1,6),'max economic depth','BackgroundColor',[1 1 1],'Margin',1)

saveas(gcf,fullfile('detailedOutput',parameterSetID{1},[parameterSetID{2} '_levelCost']),'epsc')

waterUseFig = figure();
set(waterUseFig,'PaperPosition',[1 1 3.25 6])
myTechColors = [141 211 199; 255 255 179; 190 186 218]/255;
waterTotal = water.*shares;
titles = {'CP Flat Tariff','CP Real Cost','Optimal Management'};
for ii=1:3
    subtightplot(4,1,ii,.05,.1,[.2 .1])
    thisPlot = area(waterTotal(:,:,4-ii));
    for jj=1:numel(thisPlot)
        set(thisPlot(jj),'FaceColor',myTechColors(jj,:));
    end
    axis([1 plotYears 0 50])
    ylabel('1000 m^{3}/ha')
    title(titles{ii})
    if ii<3
        set(gca,'XTickLabel',{})
    else
        set(gca,'XTick',[1 25 50],'XTickLabel',{num2str(1) num2str(25) num2str(50)})
        xlabel('Years')
        myLegend = legend('Tech 1','Tech 2','Tech 3');
        set(myLegend,'Orientation','Horizontal')
        myPos = get(myLegend,'Position');
        myPos(1) = (1-myPos(3))/2;
        myPos(2) = .25-2*myPos(4);
        set(myLegend,'Position',myPos)
    end
    
end
saveas(waterUseFig,fullfile('detailedOutput',parameterSetID{1},[parameterSetID{2} '_waterArea']),'epsc')
