function [output] = solveCase(P,modelOpts,parameterSetID)

%% compute derived parameter values
P.h0 = P.landHeight - P.initialLift;

P.investInd = 1; P.gwDugInd = 2; P.gwBoreInd = 3; P.sbInd = 1; P.levelInd = 2; P.yrsInd = 3;

%per unit cost function for dug wells has form C(h) = ae^(bh). Set so its
%value is equal to the demand intercept at maxDepthDug and its slope is
%equal to slopeMaxDepth at this same value.
P.costDug_a = P.dDugInt*exp(-P.maxDepthDug*P.slopeMaxDepth/P.dDugInt); P.costDug_b = P.slopeMaxDepth/P.dDugInt;

%add a cost as we near the bottom of the aquifer that ensures it is
%uneconomic to pump water at the bottom of the aquifer.
P.costBore_a = log(P.dBoreInt);
P.costBore_b = P.slopeBoreZero*exp(-P.costBore_a);
P.costBore_a = 0;

%scaling parameters
%P.valScale = 1/2*max(P.dDugInt^2/P.dDugSlope,P.dBoreInt^2/P.dBoreSlope);

%% solve the optimal management problem
model.func = 'optFunc';
model.discount = P.discount;
model.params = {P};

smin(P.levelInd) = P.bottom;
%smin(P.levelInd) = 10;
smax(P.levelInd) = P.landHeight;
smin(P.sbInd) = P.shrBore0;
smax(P.sbInd) = 1;

n = [modelOpts.capNodes modelOpts.heightNodes];
fspace = fundefn('lin',n,smin,smax);
%make sure there's a break point at max depth dug
% heightBreaks = fspace.parms{P.levelInd}{1};
% if ~any(heightBreaks==P.maxDepthDug)
%     fspace.parms{P.levelInd}{1} = sort([heightBreaks; P.maxDepthDug]);
% end

snodes = funnode(fspace);
s = gridmake(snodes);

%make guess at the optimal actions
lifts = P.landHeight - s(:,P.levelInd);
costBore = P.electricity*lifts;
costDug = P.costDug_a*exp(P.costDug_b*lifts);

states = length(s);
[lb,ub] = feval(model.func,'b',s,ones(states,3),[],P);

xGuess(:,P.gwDugInd) = max((P.dDugInt - costDug)./P.dDugSlope,0);
xGuess(:,P.gwBoreInd) = max((P.dBoreInt - costBore)./P.dBoreSlope,0);
nb = netBen(xGuess(:,P.gwDugInd),xGuess(:,P.gwBoreInd),xGuess(:,P.investInd),s(:,P.levelInd),s(:,P.sbInd),P);
deltaBen = nb.bore - nb.dug;
newShr = max(s(:,P.sbInd),normcdf(deltaBen*P.discount/(1-P.discount),P.investCostMean,P.investCostSD));

xGuess = .9*xGuess;

vGuess = feval(model.func,'f',s,xGuess,[],P);
nextStateGuess = feval(model.func,'g',s,xGuess,[],P);
xGuess(:,P.investInd) = min(.5*ub(:,P.investInd),newShr-s(:,P.sbInd));

optset('dpsolve','algorithm','funcit');
optset('dpsolve','maxit',5000);
optset('dpsolve','showiters',1);
optset('dpsolve','tol',modelOpts.vtol);
optset('dpsolve','nres',1);

figure()
[c,scoord,v,x,resid,exf] = dpsolve(model,fspace,s,vGuess,xGuess);
if exf==0; keyboard; end
[shares,levels] = ndgrid(scoord{P.sbInd},scoord{P.levelInd});
s0(P.sbInd) = P.shrBore0;
s0(P.levelInd) = P.h0;
[ssim,xsim] = dpsimul(model,s0,modelOpts.minT,scoord,x);
close

output.opt.vFunc = v;
output.opt.val = funeval(c,fspace,s0);
output.opt.statePath = squeeze(ssim)';
output.opt.controlPath = squeeze(xsim)';
output.opt.valPath = feval(model.func,'f',output.opt.statePath,output.opt.controlPath,[],P);
output.opt.optVal = (P.discount.^(0:length(output.opt.valPath)-1))*output.opt.valPath;
%extract and store necessary optimal management output

% solve the adaptive expectations management problem
output.aeOut = aeSolve2(P,modelOpts);
output.aeOut.pgain = (output.opt.val - output.aeOut.aeVal)/output.aeOut.aeVal;

figure()
xTitles = {'Investment','Water (Traditional)','Water (Modern)'};
sTitles = {'Share in Modern Agriculture','Pumping Lift'};
sYlabel = {'%','meters'};

% solve the "rational" expectations management problem
 output.reOut = reSolve(P,modelOpts);
 output.reOut.pgain = (output.opt.val - output.reOut.reVal)/output.reOut.reVal;
 
 yrs = min([length(output.reOut.controlPath) length(output.aeOut.controlPath) length(output.opt.controlPath)]);

 for ii=2:3; 
    subplot(2,2,ii+1); 
    plot(output.aeOut.controlPath(1:yrs,ii)/10,'--'); 
    hold on; 
    plot(output.reOut.controlPath(1:yrs,ii)/10,'-.'); 
    plot(output.opt.controlPath(1:yrs,ii)/10); 
    title(xTitles{ii});
    legend('Common Property AE','Common Property RE','Optimal Management')
end; 

for ii=1:2; 
    subplot(2,2,ii); 
    if ii==2
        constant = P.landHeight; slope = -1;
    else
        constant = 0; slope=1;
    end
    plot(constant+slope*output.aeOut.statePath(1:yrs,ii),'--'); 
    hold on; 
    plot(constant+slope*output.reOut.statePath(1:yrs,ii),'-.'); 
    plot(constant+slope*output.opt.statePath(1:yrs,ii)); 
    title(sTitles{ii});
    ylabel(sYlabel{ii});
    legend('Common Property AE','Common Property RE','Optimal Management')
end;



%% compare outputs and return key results

%% store details for later reference if needed
if ~exist('detailedOutput','dir')
    mkdir('detailedOutput')
end
if ~exist(['detailedOutput/' parameterSetID.timeStamp],'dir')
    mkdir(['detailedOutput/' parameterSetID.timeStamp])
end
save(['detailedOutput/' parameterSetID.timeStamp '/' parameterSetID.case ''])
saveas(gcf,fullfile('detailedOutput',parameterSetID.timeStamp,[parameterSetID.case '_paths']),'epsc')
