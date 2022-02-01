function hog_feature_vector = hog(img)
    % Finds the HOG feature vector for any given image

    % INPUT => img (input image)
    % OUTPUT => HOG feature vector for that particular image

    % convert rgb to grayscale
    if size(img, 3) == 3
        img = rgb2gray(img);
    end

    img = double(img);
    rows = size(img, 1);
    cols = size(img, 2);
    Ix = img;                 %Basic Matrix assignment
    Iy = img;                 %Basic Matrix assignment

    % gradients in X and Y direction. Iy is the gradient in X direction and Iy is the gradient in Y direction
    for i = 1:rows-2
        Iy(i, :)=(img(i, :) - img(i+2, :));
    end
    for i = 1:cols-2
        Ix(:, i) = (img(:, i) - img(:, i+2));
    end
    gauss=fspecial('gaussian', 8);      % Initialized a gaussian filter with sigma=0.5 * block width.    
    angle=atand(Ix ./ Iy);              % Matrix containing the angles of each edge gradient
    angle=imadd(angle, 90);             % Angles in range (0, 180)
    magnitude=sqrt(Ix.^2 + Iy.^2);
    % figure, imshow(uint8(angle));
    % figure, imshow(uint8(magnitude));
    % Remove redundant pixels in an image. 
    angle(isnan(angle)) = 0;
    magnitude(isnan(magnitude)) = 0;
    hog_feature_vector = [];          %initialized the hog_feature_vector

    % Iterations for Blocks
    for i = 0:rows/8 - 2
        for j= 0:cols/8 -2

            mag_patch = magnitude(8*i+1 : 8*i+16, 8*j+1 : 8*j+16);
            %mag_patch = imfilter(mag_patch,gauss);
            ang_patch = angle(8*i+1 : 8*i+16 , 8*j+1 : 8*j+16);
            
            block_feature=[];
            
            %Iterations for cells in a block
            for x= 0:1
                for y= 0:1
                    angleA =ang_patch(8*x+1:8*x+8, 8*y+1:8*y+8);
                    magA   =mag_patch(8*x+1:8*x+8, 8*y+1:8*y+8); 
                    histr  =zeros(1, 9);
                    
                    %Iterations for pixels in one cell
                    for p=1:8
                        for q=1:8
                            alpha= angleA(p, q);
                            
                            % Binning Process (Bi-Linear Interpolation)
                            if alpha>10 && alpha<=30
                                histr(1)=histr(1)+ magA(p, q)*(30-alpha)/20;
                                histr(2)=histr(2)+ magA(p, q)*(alpha-10)/20;
                            elseif alpha>30 && alpha<=50
                                histr(2)=histr(2)+ magA(p, q)*(50-alpha)/20;                 
                                histr(3)=histr(3)+ magA(p, q)*(alpha-30)/20;
                            elseif alpha>50 && alpha<=70
                                histr(3)=histr(3)+ magA(p, q)*(70-alpha)/20;
                                histr(4)=histr(4)+ magA(p, q)*(alpha-50)/20;
                            elseif alpha>70 && alpha<=90
                                histr(4)=histr(4)+ magA(p, q)*(90-alpha)/20;
                                histr(5)=histr(5)+ magA(p, q)*(alpha-70)/20;
                            elseif alpha>90 && alpha<=110
                                histr(5)=histr(5)+ magA(p, q)*(110-alpha)/20;
                                histr(6)=histr(6)+ magA(p, q)*(alpha-90)/20;
                            elseif alpha>110 && alpha<=130
                                histr(6)=histr(6)+ magA(p, q)*(130-alpha)/20;
                                histr(7)=histr(7)+ magA(p, q)*(alpha-110)/20;
                            elseif alpha>130 && alpha<=150
                                histr(7)=histr(7)+ magA(p, q)*(150-alpha)/20;
                                histr(8)=histr(8)+ magA(p, q)*(alpha-130)/20;
                            elseif alpha>150 && alpha<=170
                                histr(8)=histr(8)+ magA(p, q)*(170-alpha)/20;
                                histr(9)=histr(9)+ magA(p, q)*(alpha-150)/20;
                            elseif alpha>=0 && alpha<=10
                                histr(1)=histr(1)+ magA(p, q)*(alpha+10)/20;
                                histr(9)=histr(9)+ magA(p, q)*(10-alpha)/20;
                            elseif alpha>170 && alpha<=180
                                histr(9)=histr(9)+ magA(p, q)*(190-alpha)/20;
                                histr(1)=histr(1)+ magA(p, q)*(alpha-170)/20;
                            end
                        end
                    end
                    block_feature = [block_feature histr]; % concatenation of Four histograms to form one block hog_feature_vector
                end
            end

            % normalize the values in the block using L1-Norm
            block_feature = block_feature/sqrt(norm(block_feature)^2+.01);

            hog_feature_vector = [hog_feature_vector block_feature]; % features concatenation
        end
    end
    hog_feature_vector(isnan(hog_feature_vector)) = 0; % removing Infinitiy values

    % normalization of the hog_feature_vector vector using L2-Norm
    hog_feature_vector=hog_feature_vector/sqrt(norm(hog_feature_vector)^2+.001);
    for z=1:length(hog_feature_vector)
        if hog_feature_vector(z)>0.2
         hog_feature_vector(z)=0.2;
     end
    end

    hog_feature_vector=hog_feature_vector/sqrt(norm(hog_feature_vector)^2+.001);

end