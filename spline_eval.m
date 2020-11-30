function [ y ] = spline_eval( x, c, d, knots )

xmin = knots(1);
xmax = knots(end);
t = [repmat(xmin, [1, d]), knots, repmat(xmax, [1, d])];

ncoeff = numel(knots) + d - 1;

B = bspline_v2( x, ncoeff, d, t);

y = B*c;
