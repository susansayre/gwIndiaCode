function pumpingLift = glLiftFit(coeffs,timeVals,lifts)

growthRate = coeffs(1);
inflectionT = coeffs(2);

h = 1./(1+exp(-r*(timeVals-inflectionT)));

error = lifts'*lifts -(h'*lifts)^2/(h'*h);


