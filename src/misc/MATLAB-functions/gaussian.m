function H = gaussian(x, y, sigma)
	% GAUSSIAN(n, sigma)
	% dimensions: 'x' and 'y' 
	% SD: 'sigma'
	sigma_sq = 2 * sigma^2;
	H = zeros(x, y);
	for ii = 1:x;
		for jj = 1:y;
			xx = ii - x/2;
			yy = jj - y/2;
			H(ii, jj) = exp(-(xx^2 + yy^2) / sigma_sq);
		end
	end
end