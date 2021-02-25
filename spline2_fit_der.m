function c = spline2_fit_der( x, y, zx, zy, d, knots_x, knots_y, lambda, c_regularized, weights)
%lambda is regularization parameter, e.g. 0.005

if nargin < 8
	lambda = 0;
end
if nargin < 9
	c_regularized = [];
end

x=x(:);
y=y(:);
zx=zx(:);
zy=zy(:);

if nargin < 10
	weights = [];
else
    weights = weights(:);
end


xmin = knots_x(1);
xmax = knots_x(end);
tx = [repmat(xmin, [1, d]), knots_x, repmat(xmax, [1, d])];

ymin = knots_y(1);
ymax = knots_y(end);
ty = [repmat(ymin, [1, d]), knots_y, repmat(ymax, [1, d])];

nc_x = numel(knots_x) + d - 1;
nc_y = numel(knots_y) + d - 1;

if numel(x)*nc_x*nc_y<0.1*1024^3/8
	Bx = zeros(numel(x), nc_x*nc_y);
	By = zeros(numel(x), nc_x*nc_y);
else
	Bx = spalloc(numel(x), nc_x*nc_y, numel(x)*(d+1)*(d+1)); % òî÷íåå ðàññ÷èàòü
	By = spalloc(numel(x), nc_x*nc_y, numel(x)*(d+1)*(d+1));
	if nnz(bspline_x)/numel(bspline_x) < 0.1
		bspline_x = sparse(bspline_x);
		bspline_der_x = sparse(bspline_der_x);
	end
	if nnz(bspline_y)/numel(bspline_y) < 0.1
		bspline_y = sparse(bspline_y);
		bspline_der_y = sparse(bspline_der_y);
	end
end

[bspline_x, bspline_der_x] = bspline_v2( x, nc_x, d, tx);
[bspline_y, bspline_der_y] = bspline_v2( y, nc_y, d, ty);

for j = 1 : nc_x
    for k = 1 : nc_y
        Bx(:, (j - 1)*nc_y + k) = bspline_der_x(:, j).*bspline_y(:, k);
        By(:, (j - 1)*nc_y + k) = bspline_x(:, j).*bspline_der_y(:, k);
    end
end
assert(nnz(Bx) <= numel(x)*(d+1)*(d+1))


if ~isempty(weights)
    BxT = Bx'*diag(sparse(weights));
    ByT = By'*diag(sparse(weights));
else
    BxT = Bx';
    ByT = By';
end

%c = (Bx'*Bx+By'*By )\(Bx'*zx + By'*zy);
if numel(x)*nc_x*nc_y<0.1*1024^3/8
    if isempty(c_regularized)
    	c = (BxT*Bx+ ByT*By + lambda*eye(nc_x*nc_y)) \ (BxT*zx + ByT*zy);
    else
        c = (BxT*Bx+ ByT*By + lambda*eye(nc_x*nc_y)) \ (BxT*zx + ByT*zy + lambda * c_regularized(:));
    end
else
    if isempty(c_regularized)
        c = (BxT*Bx+ ByT*By + lambda*speye(nc_x*nc_y)) \ (BxT*zx + ByT*zy);
    else
        c = (BxT*Bx+ ByT*By + lambda*speye(nc_x*nc_y)) \ (BxT*zx + ByT*zy + lambda * c_regularized(:));
    end
end
