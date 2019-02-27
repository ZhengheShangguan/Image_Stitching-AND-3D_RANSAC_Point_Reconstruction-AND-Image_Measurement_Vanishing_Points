%% CS 543 Assignment 03
%% Part 01
%% Zhenghe's Code
%% House Clean work
close all
clear all 
clc

%% 1. Initial work: Load the Image and Convert it into gray & double format
% name of the input file
imname01 = 'uttower_left.JPG';
imname02 = 'uttower_right.JPG';
% read in the image
fullim01 = imread(imname01);
fullim02 = imread(imname02);
% get rgb channels for the final RGB compositd result
fullim01_r = fullim01(:,:,1);
fullim01_g = fullim01(:,:,2);
fullim01_b = fullim01(:,:,3);
fullim02_r = fullim02(:,:,1);
fullim02_g = fullim02(:,:,2);
fullim02_b = fullim02(:,:,3);
% convert images to grayscale and double
fullim01 = im2double(rgb2gray(fullim01));
fullim02 = im2double(rgb2gray(fullim02));
% normalize the value to [0,1]
fullim01 = fullim01/max(max(fullim01));
fullim02 = fullim02/max(max(fullim02));
% get the size of the image
[length01, width01] = size(fullim01);
[length02, width02] = size(fullim02);

%% 2. Detect feature points with Harris detector
sigma = 2;
thresh = 5e-2;
radius = 5;
disp = 0;
[cim01, row01, col01] = harris(fullim01, sigma, thresh, radius, disp);
[cim02, row02, col02] = harris(fullim02, sigma, thresh, radius, disp);

%% 3. Extract local neighborhoods of keypoints & form descriptors by flattening
% neighborhood size
neigh = 7;
% eliminate the points too close to the edges
[num01,~] = size(row01);
[num02,~] = size(row02);
list01 = [];
list02 = [];
for i = 1:num01
    if row01(i) < (neigh+1)/2 | row01(i) > length01-(neigh-1)/2 | col01(i) < (neigh+1)/2 | col01(i) > width01-(neigh-1)/2
        list01 = [list01,i];
    end
end
for i = 1:num02
    if row02(i) < (neigh+1)/2 | row02(i) > length02-(neigh-1)/2 | col02(i) < (neigh+1)/2 | col02(i) > width02-(neigh-1)/2
        list02 = [list02,i];
    end
end        
row01(list01,:) = [];
col01(list01,:) = [];
row02(list02,:) = [];
col02(list02,:) = [];
[num01,~] = size(row01);
[num02,~] = size(row02);
    
% initialize the descriptor array
descr01 = zeros(num01, neigh^2);
descr02 = zeros(num02, neigh^2);
% calculate the original descriptors
for i = 1:num01
    for j = 1:neigh^2
        [m,n] = ind2sub([neigh,neigh],j);
        descr01(i,j) = fullim01(row01(i)+m-(neigh+1)/2,col01(i)+n-(neigh+1)/2);
    end
end
for i = 1:num02
    for j = 1:neigh^2
        [m,n] = ind2sub([neigh,neigh],j);
        descr02(i,j) = fullim02(row02(i)+m-(neigh+1)/2,col02(i)+n-(neigh+1)/2);
    end
end

%% 4. Compute distances btw descriptors
% normalize the descriptors first
for i = 1:num01
    descr01(i,:) = descr01(i,:) - mean(descr01(i,:));
    descr01(i,:) = descr01(i,:) / norm(descr01(i,:));
end
for i = 1:num02
    descr02(i,:) = descr02(i,:) - mean(descr02(i,:));
    descr02(i,:) = descr02(i,:) / norm(descr02(i,:));
end
% compute the Euclidean distances
dist_eu = dist2(descr01, descr02);

%% 5. Select putative all matches below a threshold
thre_dist = 1e-1;
dist_refined = dist_eu .* (dist_eu<thre_dist);
[row,col] = find(dist_refined);

%% 6. Run RANSAC to estimate a homography
% correspondingly put the original putative matches' position into a single matrix
[num_dsc,~] = size(row);
matches01 = zeros(num_dsc,2);
matches02 = zeros(num_dsc,2);
% for an image coordinates: x-axis represents column, y-axis represents row.
for i = 1:num_dsc
    matches01(i,1) = col01(row(i)); 
    matches01(i,2) = row01(row(i));
end
for i = 1:num_dsc
    matches02(i,1) = col02(col(i));
    matches02(i,2) = row02(col(i));
end
% RANSAC parameter configuration
iterNum = 20000;
thre_rsc = 3;
sampleNum = 4;
% RANSAC optimal result initialization
iterOpt = 0;
inlierNumOpt = 0;
avrResOpt = 0;
H_Opt = zeros(3,3);
inliers = cell(1,iterNum);
% run the RANSAC algorithm
for i = 1:iterNum
    % Sampling matches (randomly choose 4)
    sampleId = randperm(num_dsc,sampleNum);
    pts1 = matches01(sampleId,:);
    pts2 = matches02(sampleId,:);
    
    % find an initial fitting homography based on the samples
    H = findHomo(pts1, pts2, sampleNum);
    
    % run my RANSAC function
    [inlierNum, avrRes, inliers{i}] = my_ransac(H, matches01, matches02, thre_rsc);
    
    if inlierNumOpt < inlierNum
        inlierNumOpt = inlierNum;
        avrResOpt = avrRes;
        iterOpt = i;
        H_Opt = H;
    end
end
inlierOpt = inliers{iterOpt};

% %% temperory part: show the matching results for final inliers
% imshow([fullim01 fullim02]); hold on;
% plot(inlierOpt(:,1), inlierOpt(:,2), '+r');
% plot(inlierOpt(:,3)+size(fullim01,2), inlierOpt(:,4), '+r');
% line([inlierOpt(:,1) inlierOpt(:,3) + size(fullim01,2)]', inlierOpt(:,[2 4])', 'Color', 'r');
% pause;

%% report the number, the average residual, inlier locations on the plot
figure, imagesc(fullim01), axis image, colormap(gray), hold on
plot(inlierOpt(:,1),inlierOpt(:,2),'bs'), title('RANSAC inlier matches-left image'),
legend(['inlier number = ' num2str(inlierNumOpt) ', average Residuals = ' num2str(avrResOpt)]);
figure, imagesc(fullim02), axis image, colormap(gray), hold on
plot(inlierOpt(:,3),inlierOpt(:,4),'bs'), title('RANSAC inlier matches-right image');
legend(['inlier number = ' num2str(inlierNumOpt) ', average Residuals = ' num2str(avrResOpt)]);

%% 7. Warp one image onto the other with the estimated transformation
% create a homography transformation structure, we need a transpose for
% Homography because for the image: x-axis represents column, y-axis
% represents row.
Tform = maketform('projective',H_Opt');
% apply the transformation to the left image, I need the translation info
[im_new01,xdata, ydata] = imtransform(fullim01,Tform,'bicubic');
% figure, imshow(im_new01,'XData',xdata,'YData',ydata), axis on, axis([xdata, ydata]), colormap(gray);

%% 8. Create a new image to hold the panorama and composite the two images into it
% first, right image is needed to be translated with left_image translation
% information when give the transformation
im_new02 = imtranslate(fullim02,[-xdata(1), -ydata(1)],'OutputView','full');
% extend the matrix im_new01 with zeros
[row_final, col_final] = size(im_new02);
[row_im1, col_im1] = size(im_new01);
im_new01 = [im_new01, zeros(row_im1,col_final-col_im1);zeros(row_final-row_im1,col_final)];
% get a final grey-scale composited image
C = imfuse(im_new01,im_new02,'blend','Scaling','joint');
figure, imagesc(C), axis off, colormap(gray);

% finally, get a RGB composited image
[im_new01_R,xdata, ydata] = imtransform(fullim01_r,Tform,'bicubic');
[im_new01_G,xdata, ydata] = imtransform(fullim01_g,Tform,'bicubic');
[im_new01_B,xdata, ydata] = imtransform(fullim01_b,Tform,'bicubic');

im_new02_R = imtranslate(fullim02_r,[-xdata(1), -ydata(1)],'OutputView','full');
im_new02_G = imtranslate(fullim02_g,[-xdata(1), -ydata(1)],'OutputView','full');
im_new02_B = imtranslate(fullim02_b,[-xdata(1), -ydata(1)],'OutputView','full');

im_new01_R = [im_new01_R, zeros(row_im1,col_final-col_im1);zeros(row_final-row_im1,col_final)];
im_new01_G = [im_new01_G, zeros(row_im1,col_final-col_im1);zeros(row_final-row_im1,col_final)];
im_new01_B = [im_new01_B, zeros(row_im1,col_final-col_im1);zeros(row_final-row_im1,col_final)];
C_R = imfuse(im_new01_R,im_new02_R,'blend','Scaling','joint');
C_G = imfuse(im_new01_G,im_new02_G,'blend','Scaling','joint');
C_B = imfuse(im_new01_B,im_new02_B,'blend','Scaling','joint');
C_RGB = cat(3,C_R,C_G,C_B);

figure, imshow(C_RGB), axis off;



