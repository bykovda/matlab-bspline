function [bspline_x, bspline_y, B] = spline2_basis( x, y, d, knots_x, knots_y)

if nargin < 7
	lambda = 0;
end

x = x(:);
y = y(:);

xmin = knots_x(1);
xmax = knots_x(end);
tx = [repmat(xmin, [1, d]), knots_x, repmat(xmax, [1, d])];

ymin = knots_y(1);
ymax = knots_y(end);
ty = [repmat(ymin, [1, d]), knots_y, repmat(ymax, [1, d])];

ncoeff_x = numel(knots_x) + d - 1;
ncoeff_y = numel(knots_y) + d - 1;

%{
bspline_x = zeros(numel(x), ncoeff_x);
bspline_y = zeros(numel(x), ncoeff_y);
parfor j = 1 : ncoeff_x
	bspline_x(:, j) = bspline(x, j, d, tx, d);
end
parfor k = 1 : ncoeff_y
	bspline_y(:, k) = bspline(y, k, d, ty, d);
end
%}

bspline_x = bspline_v2( x, ncoeff_x, d, tx);
bspline_y = bspline_v2( y, ncoeff_y, d, ty);

if nargout == 3
	B = zeros(numel(x), ncoeff_x*ncoeff_y);
	for j = 1 : ncoeff_x
		for k = 1 : ncoeff_y
			B(:, (j - 1)*ncoeff_y + k) = bspline_x(:, j).*bspline_y(:, k);
		end
	end
end
	


