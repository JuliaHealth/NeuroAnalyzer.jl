function [m_norm, z_norm, pp1, pp2] = compare_images(img1, img2)
	% COMPARE_IMAGES(IMG1, IMG2)
	% Calculate the difference and its norms.
	% Output:
	% 'm_norm' Manhattan norm
	% 'z_norm' Zero norm: how many non-zero elements
	% pp1: img1 size
	% pp2: img2 size

	difference = img1 - img2;
	m_norm = sum(abs(difference), 'all');
	z_norm = sum(difference ~= 0, 'all');
	[x, y] = size(img1);
	pp1 = x * y;
	[x, y] = size(img2);
	pp2 = x * y;
end