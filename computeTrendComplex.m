function [yrsLeft,params_tau,hitsBottomAt] = computeTrend(levelSeries,levelYears);
    %takes a series of consecutive gw levels and computes the years from
    %t=0 until a steady-state is reached and the parameters for the change as a
    %function of years remaining
	order = length(levelSeries)-1;
    if ~exist('levelYears')
        levelYears = 0:order;
    end
    params_t = polyfit(levelYears,levelSeries,order);
    dHParams_t = polyder(params_t);
    Tvals = roots(dHParams_t);
    realTvals = Tvals(find(imag(Tvals)==0));
    posTvals = realTvals(find(realTvals>0));
    steadyStates = polyval(params_t,posTvals);
    internalSS = find((steadyStates>0).*(steadyStates<1));
    
    if any(internalSS)
        yrsLeft = min(posTvals(internalSS));
        tau = yrsLeft - levelYears;
        params_tau = polyfit(tau,levelSeries,order);
        hitsBottomAt = [];
    else
        disp('No internal steady state')
        %find the point where the aquifer hits the bottom
        bottomTvals=roots(params_t);
        realRoots = bottomTvals(find(imag(bottomTvals)==0));
        yrsLeft = min(realRoots(find(realRoots>0)));
        
        if yrsLeft
            hitsBottomAt = yrsLeft;
            tau = yrsLeft - levelYears;
            params_tau = polyfit(tau,levelSeries,order);
        else
            disp('cubic doesn''t work')
            
    end