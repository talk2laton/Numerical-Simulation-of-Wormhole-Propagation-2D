function xm = harmmean(X, dim)
if nargin < 2
    dim = 1;
end
xm = (mean(X.^-1, dim).^-1);