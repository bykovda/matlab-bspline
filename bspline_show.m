function bspline_show( x, d, knots )

xmin = knots(1);
xmax = knots(end);
t = [repmat(xmin, [1, d]), knots, repmat(xmax, [1, d])];

[ f, df ] = bspline_v2( x(:), ncoeff, d, t );

figure
hold all
for j = 1 : ncoeff
    plot(x(:), f(:,j));
end

figure
hold all
for j = 1 : ncoeff
    plot(x(:), df(:,j));
end
