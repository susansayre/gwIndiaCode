plotYears = 50;
manageScens = {'reCp' 'reRealE' 'opt'};

dFs = P.discount.^(0:1:plotYears-1)';
for ii=1:3
eval(['investCosts(:,ii) = dFs.*output.' manageScens{ii} '.nbDetailPath(1:plotYears,4);'])
eval(['tariffs(:,ii) = dFs.*output.' manageScens{ii} '.nbDetailPath(1:plotYears,7);'])
eval(['energyCost(:,ii) = dFs.*output.' manageScens{ii} '.nbDetailPath(1:plotYears,2);'])
eval(['benefits(:,ii) = dFs.*output.' manageScens{ii} '.nbDetailPath(1:plotYears,1);'])
end

costArray(:,4,:) = investCosts;
costArray(:,3,:) = energyCost;
costArray(:,3,1) = tariffs(:,1); 
costArray(:,2,1) = energyCost(:,1) - tariffs(:,1); %subsidy under reCp
costArray(:,1,:) = benefits - investCosts - energyCost;

accumulatedCost = cumsum(costArray,1);
costFig = figure();
set(costFig,'PaperPosition',[1 1 3.25 6],'PaperPositionMode','auto')
titles = {'CP Flat Tariff','CP Real Cost','Optimal Management'};
myColors = [43 131 186; 253 174 97; 253 174 97; 215 25 28]/255;
hatchStyle = {'fill' 'single' 'fill' 'fill'};
for ii=1:3
    subtightplot(4,1,ii,.05,.1,[.2 .1])
    thisPlot = area(accumulatedCost(:,:,ii));
    for jj=1:numel(thisPlot)
        set(thisPlot(jj),'FaceColor',myColors(jj,:));
     end
    axis([1 plotYears 0 750])
    ylabel('k Rs./ha')
    set(gca,'YTick',[0 250 500 750])
    title(titles{ii})
    if ii<3
        set(gca,'XTickLabel',{})
    else
        set(gca,'XTick',[1 25 50],'XTickLabel',{num2str(1) num2str(25) num2str(50)})
        xlabel('Years')
        [myLegend1,myLegend2,~,~] = legendflex(thisPlot,{'Social Net Benefit','Energy Cost Subsidy','Energy Cost','Investment and Maintenance Costs'},'anchor',[6 2],'buffer',[0 -25]);
        for jj=1:4
            hatchfill2(myLegend2(4+jj).Children,hatchStyle{jj})
        end
    end
    for jj=1:numel(thisPlot)
        hatchfill2(thisPlot(jj),hatchStyle{jj})
    end
end
saveas(costFig,fullfile('detailedOutput',parameterSetID{1},[parameterSetID{2} '_costArea']),'epsc')