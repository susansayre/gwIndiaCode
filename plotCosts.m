
areaOrder = {'profit' 'waterCost' 'tariff'};
lineOrder = {'revenue' 'cost'};
caseOrder = {'opt' 'reRealE' 'reCp' 'reCp'};

for ii=1:numel(caseOrder)
    if ii==numel(caseOrder)
       farmCase = 'farmDetailI';
    else
       farmCase = 'farmDetailS';
    end
    for jj=1:numel(areaOrder)
        eval(['plotArray(:,jj,:,ii) = output.' caseOrder{ii} '.' farmCase '.' areaOrder{jj} '.*(output.' caseOrder{ii} '.sharePath>1e-3);'])
    end
    for jj=1:numel(lineOrder)
        switch lineOrder{jj}
            case 'revenue'
                eval(['lineArray(:,jj,:,ii) = output.' caseOrder{ii} '.' farmCase '.rev.*(output.' caseOrder{ii} '.sharePath>1e-3);'])
            case 'cost'
                eval(['lineArray(:,jj,:,ii) = (output.' caseOrder{ii} '.' farmCase '.waterCost + output.' caseOrder{ii} '.' farmCase '.tariff).*(output.' caseOrder{ii} '.sharePath>1e-3);'])
        end
    end
end

caseNames = {'opt' 'reRealE' 'reCp-Social' 'reCp-Farm'};
figure()
for ii=1:numel(caseNames)
    for jj=1:P.numTech
        subplot(P.numTech,numel(caseNames),(jj-1)*numel(caseNames)+ii)
        area(plotArray(:,:,jj,ii))
        title([caseNames{ii} ' Tech ' num2str(jj)])
    end
end

figure()
for ii=1:numel(caseNames)
    for jj=1:P.numTech
        subplot(P.numTech,numel(caseNames),(jj-1)*numel(caseNames)+ii)
        plot(lineArray(:,:,jj,ii))
        title([caseNames{ii} ' Tech ' num2str(jj)])
    end
end

