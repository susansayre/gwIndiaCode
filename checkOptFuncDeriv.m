%check opt func derivatives at current values

[ns,ds] = size(s); 
[nsCheck,dx] = size(x);

if ns-nsCheck
    keyboard;
end

[baseV,base_dV,base_ddV] = optFunc('f',s,x,[],P);

deltaX = 1e-4;
for ii=1:dx;
    xhat = x;
    xhat(:,ii) = xhat(:,ii) + deltaX;
    [newV,new_dV] = optFunc('f',s,x,[],P);
    deltaV(:,ii) = newV - baseV;
    delta_dV(:,:,ii) = new_dV - base_dV;
end

dV_compute = deltaV/deltaX;
ddV_compute = delta_dV/deltaX;