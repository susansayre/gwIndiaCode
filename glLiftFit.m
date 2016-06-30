function [error,K] = glLiftFit(coeffs,timeVals,lifts,nu)

growthRate = coeffs(1);
t0 = coeffs(2);

if numel(coeffs)==3
    if numel(timeVals)<3
        error('Can''t fit 3 coefficients with only two data points')
    else
        nu = coeffs(3);
    end
else
    if ~exist('nu','var')
        nu = 1;
    end
end

h = 1./((1+exp(-growthRate*nu*(timeVals-t0))).^(1/nu));

error = lifts'*lifts -(h'*lifts)^2/(h'*h);
K = (h'*lifts)/(h'*h);
%keyboard


