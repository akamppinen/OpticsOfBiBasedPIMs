function I = numInt(x1,x2,xData,yData)
%numInt performs numerical integration with Matlab's trapz function. It
%pads inner points x with starting (x1) and ending (x2) points and linearly
%extrapolates y data from inner points to the edges.

% Check inputs
if size(xData,1)~=1 || size(xData,2) < 2
    error('Give x as a row vector')
elseif size(xData)~=size(yData)
    error('x and y data needs to be of same size.')
end

% Sort data
[xData,sortInds] = sort(xData,'ascend');
yData = yData(sortInds);

% Calculation
if x1<x2
    xRange = xData>x1 & xData<x2;
    X = [x1,xData(xRange),x2];
elseif x2<x1
    xRange = xData>x2 & xData<x1;
    X = [x2,xData(xRange),x1];
    yData = -yData;
end
Y = interp1(xData(xRange),yData(xRange),X,'linear','extrap');
I = trapz(X,Y);

end
