function blp = butter_lp(img, d, n)
	% BUTTER_HP(img, d, n)
	% LP Butterworth filter of image img, cut-off d and order n.
  blp = 1 - butter_hp(img, d, n);
end