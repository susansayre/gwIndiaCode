for ci=1:12
    load(['detailedOutput/' runID '/case' num2str(ci) 'niP0'])
    output.reOut = reSolve(P,modelOpts,P.h0*ones(modelOpts.minT,1));
    output.reOut.pgain = (output.opt.val - output.reOut.reVal)/output.reOut.reVal;

    %figure('Visible',modelOpts.figureVisible)
    xTitles = {'Investment','(Traditional)','(Modern)'};
    sTitles = {'Share in Modern Agriculture','Pumping Lift'};
    sYlabel = {'%','meters'};
    % 
    yrs = min([length(output.reOut.controlPath) length(output.opt.controlPath)]);

     byFarmFig = figure();
     for ii=2:3; 
        subplot(2,2,ii-1); 
        hold on; 
        plot(output.reOut.controlPath(1:yrs,ii),'-.'); 
        plot(output.opt.controlPath(1:yrs,ii)); 
        title(['Water ' xTitles{ii}]);
        subplot(2,2,ii+1);
        hold on;
        plot(output.reOut.valPath(1:yrs,ii-1),'-.'); 
        plot(output.opt.valPath(1:yrs,ii-1));
        title(['Current Payoff ' xTitles{ii}])

    end; 

    totFig = figure();
    for ii=1:2; 
        subplot(2,2,ii); 
        if ii==2
            constant = P.landHeight; slope = -1;
        else
            constant = 0; slope=1;
        end
        hold on; 
        plot(constant+slope*output.reOut.statePath(1:yrs,ii),'-.'); 
        plot(constant+slope*output.opt.statePath(1:yrs,ii)); 
        title(sTitles{ii});
        ylabel(sYlabel{ii});
    end;

    subplot(2,2,3); 
    hold on; 
    plot(output.reOut.statePath(1:yrs,1).*output.reOut.controlPath(1:yrs,3)+(1-output.reOut.statePath(1:yrs,1)).*output.reOut.controlPath(1:yrs,2),'-.'); 
    plot(output.opt.statePath(1:yrs,1).*output.opt.controlPath(1:yrs,3)+(1-output.opt.statePath(1:yrs,1)).*output.opt.controlPath(1:yrs,2)); 
    title('Total Water Use');

    subplot(2,2,4); 
    hold on; 
    plot(output.reOut.valPath(1:yrs,3),'-.'); 
    plot(output.opt.valPath(1:yrs,3)); 
    title('Total Current Value');

    %% store details for later reference if needed
    saveas(totFig,fullfile('detailedOutput',parameterSetID{1},[parameterSetID{2} '_AggPaths']),'epsc')
    saveas(byFarmFig,fullfile('detailedOutput',parameterSetID{1},[parameterSetID{2} '_FarmPaths']),'epsc')
    clear totFig byFarmFig
    save(['detailedOutput/' parameterSetID{1} '/' parameterSetID{2} ''])
    updatedResults{ci} = output;
end

resultsniP0 = updatedResults;
save(['detailedOutput/' parameterSetID{1} '/fullResults'],'resultsniP0','-append')