function [ f, df ] = bspline_v2( x, jmax, d, t )
	cur_f = zeros(numel(x), jmax+d);
	
    %cur_d == 0
	for j = 1:jmax+d
		if (j + 1) >= numel(t) - d
			cur_f(x >= t(j), j) = 1;
		else
			cur_f(x >= t(j) & x < t(j + 1), j) = 1;
		end
	end
	
	for cur_d = 1:d
		last_f = cur_f;
		cur_f = zeros(numel(x), jmax+d);
		for j = 1:jmax+d-cur_d
			if t(j + cur_d) ~= t(j)
				leftc = (x - t(j)) / (t(j + cur_d) - t(j));
				cur_f(:, j) = cur_f(:, j) + leftc .* last_f(:, j);
			end
			if t(j + cur_d + 1) ~= t(j + 1)
				rightc = (t(j + cur_d + 1) - x) / (t(j + cur_d + 1) - t(j + 1));
				cur_f(:, j) = cur_f(:, j) + rightc .* last_f(:, j+1);
			end		
		end
	end
	
	f = cur_f(:, 1:jmax);
	if nargout == 2
		clear cur_f
		df = zeros(numel(x), jmax);
		for j = 1:jmax
			if t(j + d +1) == t(j)
				df(:,j) = zeros(size(x));
			elseif t(j) < t(j + d) && t(j + 1) == t(j + 1 +d)    
				df(:,j) = last_f(:, j) / (t(j + d)- t(j));
			elseif t(j) == t(j + d) && t(j + 1) < t(j + 1 +d)
				df(:,j) = -last_f(:, j+1) / (t(j + 1 + d) - t(j + 1));
			else 
				tmp1 =  last_f(:, j) / (t(j + d) - t(j));
				tmp2 = -last_f(:, j+1) / (t(j + d + 1) - t(j + 1));
				df(:,j) = tmp1 + tmp2;
			end
		end
		df = d * df;
	end
end
