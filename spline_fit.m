function [ c ] = spline_fit( x, y, d, knots, lambda)
%lambda is regularization parameter e.g. 0.005

if nargin < 5
	lambda = 0;
end

xmin = knots(1);
xmax = knots(end);
t = [repmat(xmin, [1, d]), knots, repmat(xmax, [1, d])];

ncoeff = numel(knots) + d - 1;

B = bspline_v2( x, ncoeff, d, t);

if lambda ~= 0
	c = (B'*B + lambda*eye(ncoeff))\(B'*y);
else
	c = B\y; 
end
