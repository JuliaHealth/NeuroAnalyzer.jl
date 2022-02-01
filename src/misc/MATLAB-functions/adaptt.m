function img_bw = adaptt(img)
	% ADAPT(IMG)
	% adaptive thresholding

	if std2(img) < 1
		img_bw = ones(size(img, 1), size(img, 2));
	else
		img_bw = im2bw(img, graythresh(img));
	end
end