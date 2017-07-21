%pick everything that will vary in any observation and identify an order
%[yrs,current tech,next tech,sol case]
clear levels shares0 shares water activeTech rawInvest
shrThresh = 1e-3;
yrInds = (1:yrs)';
scens = {'Optimal Managment',['CP Real Cost (' num2str(output.reRealE.pgain(1)*100,'%2.0f') '% gain)'], ['CP Flat Tariff (' num2str(output.reCp.pgain(1)*100,'%2.0f') '% gain)']};

%extract variables that are naturally yrs x manageCases
levels(:,:,1) = output.opt.levelPath(yrInds,:);
levels(:,:,2) = output.reRealE.levelPath(yrInds,:);
levels(:,:,3) = output.reCp.levelPath(yrInds,:);

%extract varaibles that are naturally yrs x nTech x manageCases
shares(:,:,3) = output.reCp.sharePath(yrInds,:);
shares(:,:,2) = output.reRealE.sharePath(yrInds,:);
shares(:,:,1) = output.opt.sharePath(yrInds,:);

%extract variables that are naturally manageCases x [social farms] x 1
npv(1,:) = output.opt.npv;
npv(2,:) = output.reRealE.npv;
npv(3,:) = output.reCp.npv;

pgain(1,:) = output.reRealE.pgain;
pgain(2,:) = output.reCp.pgain;

%create "tech 4"
belowThresh1 = find(-levels>P.maxDepths(1));
share1 = shares(:,1,:);
share1(belowThresh1) = 0; %set back to zero for any "fake" tech 1 use
shares(:,1,:) = share1;
activeTech = shares>shrThresh;

water(:,:,2) = output.reRealE.controlPath(yrInds,P.ind.water);
water(:,:,3) = output.reCp.controlPath(yrInds,P.ind.water);
water(:,:,1) = output.opt.controlPath(yrInds,P.ind.water);

%look at investment which is yrs x nTech x nTech x manageCases
rawInvest(:,:,3) = output.reRealE.controlPath(yrInds,P.ind.invest);
rawInvest(:,:,2) = output.reCp.controlPath(yrInds,P.ind.invest);
rawInvest(:,:,1) = output.opt.controlPath(yrInds,P.ind.invest);

% for caseInd = 1:numel(scens)
%     for ii=1:P.numTech
%         thisInvest = rawInvest(:,(1:P.numTech-1)+(ii-1)*(P.numTech-1),caseInd);
%         thisInvest = rawShr2Share(thisInvest);
%         invest(:,ii,caseInd,:) = thisInvest;
%     end
% end
myMap = 'brewer1';
numScens = numel(scens);
%plot pumping lifts
scenArray = reshape(repArray(scens,[yrs 1],2),yrs*numScens,1);

% techFig = figure();
% %plot technology shares
techIds = reshape(repmat(1:P.numTech,[yrs 1 numScens]),yrs*P.numTech*numScens,1);
scenBig = reshape(repArray(reshape(scenArray,yrs,numScens),[1 (P.numTech) 1],[1 3]),yrs*(P.numTech)*numScens,1);
% 
% techPlot = gramm('x',repmat(yrInds,(P.numTech)*numScens,1),'y',reshape(shares,yrs*(P.numTech)*numScens,1),'color',scenBig,'subset',activeTech);
% techPlot.facet_grid(techIds,[]);
% techPlot.geom_line();
% techPlot.set_names('x','Year','y','Share','color','','row','Tech');
% techPlot.set_text_options('facet_scaling',1.1);
% techPlot.set_color_options('map',myMap);
% techPlot.axe_property('YLim',[0 1]);
% %techPlot.no_legend();
% techPlot.draw();
% brewer1Color = get(techFig,'Colormap');
% 
% saveas(gcf,fullfile('detailedOutput',parameterSetID{1},[parameterSetID{2} '_techPaths']),'epsc')

combinedFig = figure();
set(combinedFig,'PaperUnits','inches');
set(combinedFig,'PaperSize',[6.5 9],'PaperPosition',[0 0 6.5 9]);

combinedPlot(1,1) = gramm('x',repmat(yrInds,numScens,1),'y',reshape(levels,yrs*numScens,1),'color',scenArray);
combinedPlot(1,1).geom_line();
combinedPlot(1,1).set_names('x','Year','y','m below surface','color','');
combinedPlot(1,1).set_color_options('map',myMap);
combinedPlot(1,1).set_title('Water Level');
for ii=1:P.numTech
    lineArray(ii) = -P.maxDepths(ii);
end
lineArray(end+1) = -max(P.idInts)*max(P.pumpEfficiency)/(P.energyNeed*min(P.kWhCost));
combinedPlot(1,1).geom_hline('yintercept',lineArray);

waterActive = activeTech;
combinedPlot(1,2) = gramm('x',repmat(yrInds,(P.numTech)*numScens,1),'y',reshape(water,yrs*(P.numTech)*numScens,1),'color',scenBig,'subset',waterActive);
combinedPlot(1,2).facet_grid(techIds,[]);
combinedPlot(1,2).geom_line();
combinedPlot(1,2).set_names('x','Year','y','cubic meters/ha','color','','row','Tech');
combinedPlot(1,2).no_legend();
combinedPlot(1,2).set_text_options('facet_scaling',1.1);
combinedPlot(1,2).set_color_options('map',myMap);
combinedPlot(1,2).set_title('Water Use');

% combinedPlot(2,2) = gramm('x',reshape(repmat({'society' 'farmers'},numScens,1),numScens*2,1),'y',reshape(npv,numScens*2,1),'color',repmat(scens,2,1));
% combinedPlot(2,2).geom_bar('stacked',0);
% combinedPlot(2,2).no_legend;
combinedPlot.set_title(P.csString);
combinedPlot.set_text_options('big_title_scaling',.5);
combinedPlot.draw();

legHand=combinedPlot(1,1).legend_axe_handle;
currPos = get(combinedPlot(1,1).legend_axe_handle,'Position');
combinedPlot(1,1).no_legend();
levelHand = combinedPlot(1,1).facet_axes_handles(1);
currLPos = get(levelHand,'Position');
currLPos(2) = currLPos(2)+currLPos(4)/3;
currLPos(3) = .5 - currLPos(1);
currLPos(4) = currLPos(4)*2/3;

set(levelHand,'Position',currLPos);
currPos(1) = 1.1*currLPos(1);
currPos(3) = .5-currPos(1);
currPos(2) = .8*currPos(2);
currPos(4) = currLPos(2)*.8;
set(legHand,'Position',currPos);

yLims = get(levelHand,'YLim');
set(combinedFig,'currentaxes',levelHand)
for ii=1:P.numTech
    if lineArray(ii)>yLims(1) && lineArray(ii)<yLims(2)
        text(45,lineArray(ii),['Max depth for tech ' num2str(ii)])
    end
end
if lineArray(end)>yLims(1) && lineArray(ii)<yLims(2)
    text(45,lineArray(end),'Max economic depth','BackgroundColor',[1 1 1])
end
saveas(gcf,fullfile('detailedOutput',parameterSetID{1},[parameterSetID{2} '_comboFig']),'epsc')

saveas(gcf,fullfile('detailedOutput',parameterSetID{1},[parameterSetID{2} '_levelPaths']),'epsc')

areaFig = figure();
set(areaFig,'PaperUnits','inches');
set(areaFig,'PaperSize',[6.5 9],'PaperPosition',[0 0 6.5 9]);
areaPlot = gramm('x',repmat(yrInds,(P.numTech)*numScens*2,1),'y',[reshape(shares,yrs*(P.numTech)*numScens,1); reshape(water.*shares,yrs*(P.numTech)*numScens,1)],'color',[techIds; techIds]);
areaPlot.facet_grid(reshape(repmat({'Share Using Tech' 'Water Use By Tech'},yrs*(P.numTech)*numScens,1),yrs*(P.numTech)*numScens*2,1),[scenBig; scenBig],'scale','free_y');
areaPlot.geom_line();
areaPlot.set_names('x','Year','y','','color','Tech','row','','column','');
areaPlot.set_text_options('facet_scaling',.5,'big_title_scaling',.5);
areaPlot.set_color_options('map','brewer2');
areaPlot.set_order_options('column',-1);
areaPlot.axe_property('XLim',[yrInds(1) yrInds(end)]);
areaPlot.set_title(P.csString);
areaPlot.draw();

areaKids = get(areaFig,'Children');
colorOrder = get(areaKids(end),'colorOrder');
colorOrderIndex = get(areaKids(end),'colorOrderIndex');
%keyboard
legKids = get(areaPlot.legend_axe_handle,'Children');
for ii=1:(P.numTech)
    colorOrder(ii,:) = get(legKids(2*ii),'Color');
end

for ii=1:numScens
    set(areaFig,'currentaxes',areaPlot.facet_axes_handles(1,ii))
    %cla
    thisArea = area(yrInds,shares(:,:,ii));
    for jj=1:numel(thisArea)
        set(thisArea(jj),'FaceColor',colorOrder(numel(thisArea)+1-jj,:));
    end
    axis([yrInds(1) yrInds(end) 0 1])
end

for ii=1:numScens
    set(areaFig,'currentaxes',areaKids(end+1-numScens-ii))
    %cla
    thisArea = area(yrInds,shares(:,:,ii).*water(:,:,ii));
    for jj=1:numel(thisArea)
        set(thisArea(jj),'FaceColor',colorOrder(numel(thisArea)+1-jj,:));
    end
    axis([yrInds(1) yrInds(end) 0 P.maxWater])
end

saveas(gcf,fullfile('detailedOutput',parameterSetID{1},[parameterSetID{2} '_areaPaths']),'epsc')

% investFig = figure();
% % plot shares moving between technologies in the lower left
% thisSize = yrs*P.numTech^2*numScens;
% investPlot = gramm('x',repmat(yrInds,P.numTech^2*numScens,1),'y',reshape(invest,yrs*P.numTech^2*numScens,1),'lineStyle',repmat(scenBig,P.numTech,1),'subset',repmat(activeTech,P.numTech,1));
% investPlot.facet_grid(repmat(techIds,P.numTech,1),reshape(repArray(1:P.numTech,[yrs P.numTech 2 1],4),thisSize,1));
% investPlot.geom_line();
% investPlot.set_names('x','Year','y','Share Moving To','linestyle','','row','Current Tech','column','Next Tech');
% investPlot.no_legend();
% investPlot.draw();
% saveas(gcf,fullfile('detailedOutput',parameterSetID{1},[parameterSetID{2} '_investPaths']),'epsc')
% 
% clear *Fig
