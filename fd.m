function dval_dx = fd(s,x,e,varargin)
%computes numerical finited difference derivatives that are aware of a few features of the problem
%varargin should always include P as it's first element and the function we're computing derivatives of as it's last
%when we call fd with fd it will also include the inner function between P and the function

P = varargin{1};
func = varargin{end};
[ns,dx] = size(x);
baseVal = feval(func,s,x,e,varargin{1:end-1});
sizeVal = size(baseVal);
sizeOut = [sizeVal dx]; %output variable will have this size

if sizeVal(1)~=ns
    disp('There seems to be a dimension problem here in fd')
end

repValSize = ones(size(sizeVal));
repValSize(1) = dx;
repDxSize = sizeVal;
repDxSize(1) = 1;

%create a large matrix of states and actions
bigS = repmat(s,dx,1);

%when rawInvest shares are 1, the increases cause a share to exceed 1 which will cause problems. My goal is to identify
%these points and compute the derivatives at these locations by decreasing rather than increasing x.
rawInvest1 = (x==1); %ns x dx
rawInvest1(:,P.ind.water) = 0;
deltaX = P.deltaX*(1-2*rawInvest1); %ns x dx matrix with P.deltaX when normal derivs work and -P.deltaX when we need to go negative 
newX = repmat(x,[1 1 dx]);

for ii=1:dx
    newX(:,ii,ii) = newX(:,ii,ii)+deltaX(:,ii);
end
newX = reshape(permute(newX,[1 3 2]),ns*dx,dx);

newVal = feval(func,bigS,newX,e,varargin{1:end-1});

dvaldx = (newVal - repmat(baseVal,repValSize))./repmat(deltaX(:),repDxSize);

dval_dx = reshape(dvaldx,[ns dx sizeVal(2:end)]);
dvdxSize = size(dval_dx);
if numel(dvdxSize>2)
    %implies my base value has multiple dimension
    dval_dx = permute(dval_dx,[1 3:numel(sizeVal)+1 2]);
end

if ns~=size(dval_dx,1)
    disp('There seems to be a dimension problem here in fd')
end

if any(isnan(reshape(dval_dx,numel(dval_dx),1))) || any(isinf(reshape(dval_dx,numel(dval_dx),1)))
	disp('Returning nans or infs in fd')
end

