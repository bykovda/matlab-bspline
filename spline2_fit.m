function [ c ] = spline2_fit( x, y, z, d, knots_x, knots_y, lambda )
%lambda is regularization parameter e.g. 0.005

if nargin < 7
	lambda = 0;
end

x = x(:);
y = y(:);
z = z(:);

xmin = knots_x(1);
xmax = knots_x(end);
tx = [repmat(xmin, [1, d]), knots_x, repmat(xmax, [1, d])];

ymin = knots_y(1);
ymax = knots_y(end);
ty = [repmat(ymin, [1, d]), knots_y, repmat(ymax, [1, d])];

ncoeff_x = numel(knots_x) + d - 1;
ncoeff_y = numel(knots_y) + d - 1;

bspline_x = bspline_v2( x, ncoeff_x, d, tx);
bspline_y = bspline_v2( y, ncoeff_y, d, ty);

b_dense = numel(x)*ncoeff_x*ncoeff_y<0.1*1024^3/8;

if b_dense
	B = zeros(numel(x), ncoeff_x*ncoeff_y);
else
	B = spalloc(numel(x), ncoeff_x*ncoeff_y, numel(x)*(d+1)*(d+1));
end

for j = 1 : ncoeff_x
    for k = 1 : ncoeff_y
        B(:, (j - 1)*ncoeff_y + k) = bspline_x(:, j).*bspline_y(:, k);
    end
end

if b_dense
	if lambda ~= 0
		c = (B'*B + lambda*eye(ncoeff_x*ncoeff_y))\(B'*z);
	else
		c = B\z; 
	end
else
	c = (B'*B + lambda*speye(ncoeff_x*ncoeff_y))\(B'*z);	
end
