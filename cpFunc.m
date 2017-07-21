function [out1,out2,out3] = cpFunc(flag,s,x,e,P)

%in the common property problem, we have 
%-two continuous "states" 
%-- icM: individual investment cost multiplier
%-- level: water level in the aquifer
%-- neither of these is affected by the landowner's actions
%-one discrete state
%--current technology = 1,2, or 3
%-one discrete action with three possible choices
%-- 1: convert to tech 1
%-- 2: convert to tech 2 
%-- 3: convert to tech 3

    [ns,ds] = size(s);
    techs = zeros(ns,P.numTech);
    techs(sub2ind([ns P.numTech],(1:ns)',s(:,P.ind.tech)))= 1; %creates an ns x numTech matrix with ones telling us which tech the farm has
    lifts = P.landHeight - P.levelPath(s(:,P.ind.yrInd)); %ns x 1
    icMs = s(:,P.ind.icM);
    
    if any(e)
        error('This function can''t handle shocks')
    end
    
    ns = size(s,1);
    ds = size(s,2);
    dx = size(x,2);
    
	switch flag
	
   		case 'f'
            maxWaters = waterLimit(lifts,P,s(:,P.ind.tech));

            Intercepts(:,1) = P.idInts(s(:,P.ind.tech));
            Slopes(:,1) = P.idSlopes(s(:,P.ind.tech));
            eCosts(:,1) = P.eCosts(s(:,P.ind.tech));
            costs(:,1) = P.eCostShr*lifts.*eCosts; %ns x ntech
    
            water = min(maxWaters,max((Intercepts - costs)./Slopes,0));
            
            fixedECost(:,1) = P.fixedEcosts(s(:,P.ind.tech));
            fixedECost(maxWaters==0) = 0;
            
            nbFarm = Intercepts.*water - 0.5*Slopes.*(water.^2) - costs.*water - fixedECost;
            
            actionMat = zeros(ns,P.numActions);
            actionMat(sub2ind([ns P.numActions],(1:ns)',x)) = 1; %creates an ns x numActions matrix with ones denoting the action taken
            %return net benefits
            
            investCosts = sum((techs*P.techCostMat).*actionMat,2).*(1+s(:,P.ind.icM));
            out1 = nbFarm - investCosts;
 
        case 'g'
            out1(:,P.ind.yrInd) = min(s(:,P.ind.yrInd)+1,P.lastYear);
            out1(:,P.ind.icM) = icMs;

            actionMat = zeros(ns,P.numActions);
            actionMat(sub2ind([ns P.numActions],(1:ns)',x)) = 1; %creates an ns x numActions matrix with ones denoting the action taken
            
            out1(:,P.ind.tech) = sum((techs*P.transitionMat).*actionMat,2); 
                       
    end
	
    if any(isnan(out1))
        disp('cpFunc is returning nans')
	end