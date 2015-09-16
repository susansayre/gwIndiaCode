function [yrsLeft,params_tau] = computeTrend(levelSeries,levelYears);
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
    
    if isreal(Tvals)
        yrsLeft = ceil(min(Tvals'));
        if yrsLeft<0
            disp('It looks like you''re beyond the only possible steady state which is weird.')
            disp('I''m going to pretend this isn''t a problem and see what happens.')
            tau = yrsLeft - levelYears;
            params_tau = polyfit(tau,levelSeries,order);
       elseif yrsLeft==0
            disp('It looks like you''re at the steady state')
            if max(Tvals')>0
                disp('I''m going to use a future steady state')
                yrsLeft = maxTvals;
            else
                yrsLeft = 0; params_tau = [zeros(1, order) polyval(params_t,yrsLeft)];
            end
        else
            tau = yrsLeft - levelYears;
            params_tau = polyfit(tau,levelSeries,order);
        end
    end