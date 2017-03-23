function [y,p] = nancorrcoef(x1,x2)

% NANCORRCOEF computes the correlation between two time series
% NAN values are ignored

n1 = size(x1,1);
n2 = size(x2,1);
if (n1 == 1) x1=x1'; end
if (n2 == 1) x2=x2'; end

if (size(x1,1) ~= size(x2,1))
    error ('The lengths of Y and X must match.')
end

k = find (~isnan(x1) & ~isnan(x2));
x1 = x1(k);
x2 = x2(k);
[C,P] = corrcoef(x1,x2);
y = C(1,2);
p = P(1,2);
