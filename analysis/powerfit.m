function [a,b,outliers] = powerfit(x,y,RmOutlier,nMAD,OL1side)
% Fit y = a*x^b.
% RmOutlier - set to non-zero if removing outliers using scaled median absolute deviations (MAD)
% nMAD - outlier threshold, nMAD times MAD. Default: 3
% OL1side - if non-zero consider outliers on only one side depending on sign. Default: 0
if ~isequal(numel(x),numel(y))
    error('Number of input elements must agree.');
end
if any(x<=0)||any(y<=0)||~isreal(x)||~isreal(y)
    error('Input values must be real positive.')
end
if nargin<=2
    RmOutlier = 0;
else
    if nargin<=3||isempty(nMAD),    nMAD = 3;       end
    if nargin<=4||isempty(OL1side),	OL1side = 0;	end
end

X = log(x(:));
Y = log(y(:));

p = polyfit(X,Y,1);
a = exp(p(2));
b = p(1);

if RmOutlier
    err = Y - (p(1)*X+p(2));
    MAD = -1/(sqrt(2)*erfcinv(3/2))*median(abs(err-median(err)));
    if OL1side
        outliers = sign(OL1side)*err>nMAD*MAD;
    else
        outliers = abs(err)>nMAD*MAD;
    end
    p = polyfit(X(~outliers),Y(~outliers),1);
    a = exp(p(2));
    b = p(1);
    outliers = find(outliers);
else
    outliers = [];
end

end