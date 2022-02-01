function bhp = butter_hp(img, d, n)
	% BUTTER_HP(img, d, n)
	% HP Butterworth filter of image img, cut-off d and order n.
	h = size(img, 1);
	w = size(img, 2);
	[x, y] = meshgrid(-floor(w/2):floor(w-1)/2, -floor(h/2):floor(h-1)/2);
	bhp = 1./(1.0+(d./(x.^2 + y.^2).^0.5).^(2*n));
end