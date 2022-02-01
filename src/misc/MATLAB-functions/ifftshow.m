function ifftshow(imgX)
	img1 = ifft2(imgX);
	img2 = abs(img1);
	img3 = img2/max(img2(:));
	imshow(img3)
end