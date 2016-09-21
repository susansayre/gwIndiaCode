%create water use comparison figure
waterUse = panel();
waterUse.pack(1,3);

farmType = {'Traditional','Modern','Average'};
for ii=1:3
    waterUse(1,ii).select();
    plot(yrInds,controlPaths(:,1,ii+1,mInd),'r-','lineWidth',2)
    hold on;
    plot(yrInds,controlPaths(:,2,ii+1,mInd),'b--','lineWidth',2)
    ylabel('thousand m^3/ha')
    title(farmType{ii})
    axis([0 30 0 20])
end

legHandle = legend('Optimal Management','Common Property');
