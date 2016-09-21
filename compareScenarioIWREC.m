%compare solutions

%calculate differences in paths

compareCases = {'endogInvest' 'noInvestAvg'}; %cheat to make code work while waiting for results

P = endogInvest.P;
ds = size(endogInvest.opt.statePath,2);
dx = size(endogInvest.opt.controlPath,2);

yrs = 99;

for ii=1:numel(compareCases)
    eval(['theseResults = ' compareCases{ii} ';'])
    statePaths(:,1,:,ii) = theseResults.opt.statePath(1:yrs,1:2); %don't worry about invest cost path
    controlPaths(:,1,:,ii) = theseResults.opt.controlPath(1:yrs,:);
    statePaths(:,2,:,ii) = theseResults.reOut.statePath(1:yrs,:);
    controlPaths(:,2,:,ii) = theseResults.reOut.controlPath(1:yrs,:);
    valPaths(:,1,:,ii) = theseResults.opt.valPath(1:yrs,:);
    valPaths(:,2,:,ii) = theseResults.reOut.valPath(1:yrs,:);
    values(ii,1) = theseResults.opt.val;
    values(ii,2) = theseResults.reOut.reVal;
end

controlPaths(:,:,4,:) = statePaths(:,:,1,:).*controlPaths(:,:,3,:) + (1-statePaths(:,:,1,:)).*controlPaths(:,:,2,:);
valPaths(:,:,4,:) = statePaths(:,:,1,:).*valPaths(:,:,2,:) + (1-statePaths(:,:,1,:)).*valPaths(:,:,1,:); %pv payoff before investcost
valPaths(:,:,5,:) = valPaths(:,:,4,:) - valPaths(:,:,3,:); %invest costs
valPaths(:,:,6,:) = valPaths(:,:,2,:) - (1-P.discount)/P.discount*cumsum(valPaths(:,:,5,:),1); %modern payoff less annualized cost of all prior investments
valPaths(:,:,7,:) = valPaths(:,:,4,:) - (1-P.discount)/P.discount*cumsum(valPaths(:,:,5,:),1);

statePaths(:,3,:,:) = statePaths(:,2,:,:) - statePaths(:,1,:,:);
controlPaths(:,3,:,:) = controlPaths(:,2,:,:) - controlPaths(:,1,:,:);
valPaths(:,3,:,:) = valPaths(:,2,:,:) - valPaths(:,1,:,:);
mInd = 1;
%plot water use paths
yrInds = 0:yrs-1;

%create water use comparison figure
figure()

% farmType = {'Water (Traditional)','Water (Modern)','Water (Average)'};
% for ii=1:3
%     subplot(2,3,ii)
%     plot(yrInds,controlPaths(:,1,ii+1,mInd),'r-','lineWidth',2)
%     hold on;
%     plot(yrInds,controlPaths(:,2,ii+1,mInd),'b--','lineWidth',2)
%     ylabel('thousand m^3/ha')
%     title(farmType{ii})
%     axis([0 75 0 20])
% end

secondFigPlots(:,:,1) = statePaths(:,:,1,1);
secondFigPlots(:,:,2) = P.landHeight - statePaths(:,:,2,1);

pathType = {'Share in Modern', 'Pumping Lift', 'NPV of Benefits'};
ylabels={'%','meters'};
for ii=1:2
    subplot(2,2,ii)
    plot(yrInds,secondFigPlots(:,1,ii),'r-','lineWidth',2)
    hold on;
    plot(yrInds,secondFigPlots(:,2,ii),'b--','lineWidth',2)
    ylabel(ylabels{ii})
    title(pathType{ii})
    myAxis = axis;
    myAxis(2) = 60;
    axis(myAxis)
end

legHandle = legend('Optimal Management','Common Property');
set(legHandle,'Position',[0.41 0.4 0.18 0.08])

% thirdFigPlots = valPaths(1:30,:,[1 6 7],1);
% 
% valType = {'Traditional' 'Modern' 'Weighted Average'};
% for ii=1:3
%     subplot(3,3,6+ii)
%     plot(yrInds(1:30),thirdFigPlots(:,1,ii),'r-','lineWidth',2)
%     hold on;
%     plot(yrInds(1:30),thirdFigPlots(:,2,ii),'b--','lineWidth',2)
%     ylabel('$/ha')
%     title(valType{ii})
%     axis([0 30 -2 4])
%     plot(yrInds(1:30),0*yrInds(1:30),':','color',[.3 .3 .3])
% end
% 
saveas(gcf,'combinedFig','epsc')

figure()
comparisonShrLevels = [0:.01:.99];
convertCostLevels = norminv(comparisonShrLevels,P.investCostMean,P.investCostSD);
convertCostLevels(find(comparisonShrLevels<=P.shrBore0)) = 0; %conversion cost paid before investigation
statePathSize = size(statePaths);
for ii=1:numel(convertCostLevels)
    thisConvertYears = max(repmat(yrInds',[1 statePathSize(2) 1 statePathSize(4)]).*(statePaths(:,:,1,:)<comparisonShrLevels(ii)));
    thisConvertYears(find(thisConvertYears==yrs-1))=Inf;
    convertYears(ii,:,:) = thisConvertYears;
    thisValPath = valPaths(:,:,1,:).*(statePaths(:,:,1,:)<comparisonShrLevels(ii)) + valPaths(:,:,2,:).*(statePaths(:,:,1,:)>=comparisonShrLevels(ii));
    thisValPath(end+1,:,:) = thisValPath(end,:,:)*P.discount/(1-P.discount);
    byFarmValPath(:,:,:,ii) = thisValPath;
    NPVBenefits(ii,:,:)= sum(repmat(P.discount.^[yrInds'; yrInds(end)+1],[1 statePathSize(2) statePathSize(4)]).*byFarmValPath(:,:,:,ii))- convertCostLevels(ii)*P.discount.^(convertYears(ii,:,:));
end

figure()
subplot(2,1,1)
title('Net present value of benefits')
plot(comparisonShrLevels,NPVBenefits(:,1,1),'r-','lineWidth',2)
hold on;
plot(comparisonShrLevels,NPVBenefits(:,2,1),'b--','lineWidth',2)
xlabel('Investment Cost Percentile')
ylabel('NPV of Benefits (kRs/ha)')

subplot(2,1,2)
title('Percentage gain (loss)')
plot(comparisonShrLevels,(NPVBenefits(:,1,1) - NPVBenefits(:,2,1))./NPVBenefits(:,2,1)*100,'lineWidth',2)
hold on;
plot(comparisonShrLevels,0*comparisonShrLevels,':','color',[.3 .3 .3])
xlabel('Investment Cost Percentile')
ylabel('% gain (loss) over common property')

saveas(gcf,'NPVbyFarm','epsc')

for ii=2:2
    figure()
    plot(yrInds,statePaths(:,1,ii,1),'r-','lineWidth',2)
    hold on
    plot(yrInds,statePaths(:,1,ii,2),'r-.','lineWidth',1.5)
    plot(yrInds,statePaths(:,2,ii,1),'b--','lineWidth',2)
    plot(yrInds,statePaths(:,2,ii,2),'b-.','lineWidth',1.5)
    title('Pumping Lifts')
    myAxis = axis;
    myAxis(2) = 10;
    axis(myAxis)
end

figure()
plot(yrInds,controlPaths(:,1,4,1),'r-','lineWidth',2)
hold on;
plot(yrInds,controlPaths(:,1,4,2),'r-.','lineWidth',1.5)
plot(yrInds,controlPaths(:,2,4,1),'b--','lineWidth',2)
plot(yrInds,controlPaths(:,2,4,2),'b-.','lineWidth',1.5)
myAxis = axis;
myAxis(2) = 10;
axis(myAxis)
title('Total Water Use')

