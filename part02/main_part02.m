%%
%% load images and match files for the first example
%%
clc
clear all

I1 = imread('library1.jpg');
I2 = imread('library2.jpg');
P1 = load('library1_camera.txt');
P2 = load('library2_camera.txt');
matches = load('library_matches.txt'); 
% this is a N x 4 file where the first two numbers of each row
% are coordinates of corners in the first image and the last two
% are coordinates of corresponding corners in the second image: 
% matches(i,1:2) is a point in the first image
% matches(i,3:4) is a corresponding point in the second image

N = size(matches,1);

%%
%% display two images side-by-side with matches
%% this code is to help you visualize the matches, you don't need
%% to use it to produce the results for the assignment
%% 
imshow([I1 I2]); hold on;
plot(matches(:,1), matches(:,2), '+r');
plot(matches(:,3)+size(I1,2), matches(:,4), '+r');
line([matches(:,1) matches(:,3) + size(I1,2)]', matches(:,[2 4])', 'Color', 'r');
% pause;

%%
%% display second image with epipolar lines reprojected 
%% from the first image
%%

%% Choose a method for finding fundamental matrix
% method = 'unnormalized' : unnormalized algorithm to find F, also for Triangluation;
% method = 'normalized' : normalized algorithm to find F;
% method = 'RANSAC' : use ransac to find matches then find F.
method = 'normalized';


%% Normalized or Unnormalized fit_fundamental matrix

if strcmp(method,'normalized') || strcmp(method,'unnormalized')
    
    close all
    % first, fit fundamental matrix to the matches
    F = fit_fundamental(matches, method); % this is a function that you should write
    L = (F * [matches(:,1:2) ones(N,1)]')'; % transform points from 
    % the first image to get epipolar lines in the second image

    % find points on epipolar lines L closest to matches(:,3:4)
    L = L ./ repmat(sqrt(L(:,1).^2 + L(:,2).^2), 1, 3); % rescale the line
    pt_line_dist = sum(L .* [matches(:,3:4) ones(N,1)],2);
    closest_pt = matches(:,3:4) - L(:,1:2) .* repmat(pt_line_dist, 1, 2);

    % find endpoints of segment on epipolar line (for display purposes)
    pt1 = closest_pt - [L(:,2) -L(:,1)] * 10; % offset from the closest point is 10 pixels
    pt2 = closest_pt + [L(:,2) -L(:,1)] * 10;

    % display points and segments of corresponding epipolar lines
    clf;
    imshow(I2); hold on;
    plot(matches(:,3), matches(:,4), '+r');
    line([matches(:,3) closest_pt(:,1)]', [matches(:,4) closest_pt(:,2)]', 'Color', 'r');
    line([pt1(:,1) pt2(:,1)]', [pt1(:,2) pt2(:,2)]', 'Color', 'g');

    % use legend to report the mean squared distance btw points in both images 
    % and the corresponding epipolar lines
    dist_avr = mean(sum((closest_pt - matches(:,3:4)).^2,2));
    legend(['average residual = ' num2str(dist_avr)]);
    
    
    %% Particularly for Triangulation and Reconstruction (Step 4)
    % (Since normalized is unnecessary, only apply the "unnormalized" method)
    
    % find the camera center
    [U1,S1,V1] = svd(P1);
    [U2,S2,V2] = svd(P2);
    center01 = V1(:,end)/V1(end);
    center02 = V2(:,end)/V2(end);
    
    % Linear least squares to triangulate the position of the matches
    X = zeros(N,3);
    Res01_proj = zeros(N,1);
    Res02_proj = zeros(N,1);
    for i = 1:N
        skew_matrix1 = [0 -1 matches(i,2);1 0 -matches(i,1);-matches(i,2) matches(i,1) 0];
        skew_matrix2 = [0 -1 matches(i,4);1 0 -matches(i,3);-matches(i,4) matches(i,3) 0];
        Tr_tmp = [skew_matrix1 * P1; skew_matrix2 * P2];
        [U, S, V] = svd(Tr_tmp);
        tmp = (V(:,end)./V(end))';
        X(i,:) = tmp(:,1:3);
        X_proj1 = (P1 * tmp')';
        X_proj1 = X_proj1 ./ X_proj1(end);
        X_proj2 = (P2 * tmp')';
        X_proj2 = X_proj2 ./ X_proj2(end);
        Res01_proj(i) = sum((X_proj1(1:2)-matches(i,1:2)).^2,2);
        Res02_proj(i) = sum((X_proj2(1:2)-matches(i,3:4)).^2,2);
    end
    avr_proj_Res = sum(Res01_proj)/N;
    
    % Display the two camera centers and reconstructed points in 3D, report the residuals
    figure
    % subplot(1,2,1)
    hold on
    plot3(X(:,1),X(:,2),X(:,3),'b.');
    plot3(center01(1),center01(2),center01(3),'r*');
    plot3(center02(1),center02(2),center02(3),'g*');
    axis equal
    lim = axis;
    legend(['Projected Residual = ' num2str(avr_proj_Res)]);
    xlabel('X-axis')
    ylabel('Y-axis')
    zlabel('Z-axis')
    title('Reconstructed 3D points')

%     % particularly for creating a GIF image, the same thing for image02
%     create_gif(X, center01, avr_proj_Res);
    

%% Particularly for RANSAC implement (Step 3)
% when it's RANSAC, find the putative match generation and RANSAC w\o groundtruth.
elseif strcmp(method,'RANSAC')

    close all
    % 1. Initial work: Convert it into gray & double format
    % convert images to grayscale and double
    I1 = im2double(rgb2gray(I1));
    I2 = im2double(rgb2gray(I2));
    % normalize the value to [0,1]
    I1 = I1/max(max(I1));
    I2 = I2/max(max(I2));
    % get the size of the image
    [length01, width01] = size(I1);
    [length02, width02] = size(I2);
    
    % 2. Detect feature points with Harris detector
    % here I use a relatively loose threshold to get more features for selecting
    sigma = 2;
    thresh = 2e-2; % 1e-3; % 1e-1;
    radius = 5;
    disp = 0;
    [cim01, row01, col01] = harris(I1, sigma, thresh, radius, disp);
    [cim02, row02, col02] = harris(I2, sigma, thresh, radius, disp);

    % 3. Extract local neighborhoods of keypoints & form descriptors by 
    % flattening neighborhood size
    neigh = 11; % increase the neighboorhood area to make the match accurate
    % eliminate the points too close to the edges
    [num01,~] = size(row01);
    [num02,~] = size(row02);
    list01 = [];
    list02 = [];
    for i = 1:num01
        if row01(i) < (neigh+1)/2 || row01(i) > length01-(neigh-1)/2 || col01(i) < (neigh+1)/2 || col01(i) > width01-(neigh-1)/2
            list01 = [list01,i];
        end
    end
    for i = 1:num02
        if row02(i) < (neigh+1)/2 || row02(i) > length02-(neigh-1)/2 || col02(i) < (neigh+1)/2 || col02(i) > width02-(neigh-1)/2
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
            descr01(i,j) = I1(row01(i)+m-(neigh+1)/2,col01(i)+n-(neigh+1)/2);
        end
    end
    for i = 1:num02
        for j = 1:neigh^2
            [m,n] = ind2sub([neigh,neigh],j);
            descr02(i,j) = I2(row02(i)+m-(neigh+1)/2,col02(i)+n-(neigh+1)/2);
        end
    end

    % 4. Compute distances btw descriptors
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

    % 5. Select putative all matches below a threshold
    % more threshold of distance tighter to get more accurate matches
    thre_dist = 2e-1; % 2.0e-1;
    [row,col] = find(dist_eu<thre_dist);

    % 6. Run RANSAC to estimate a homography
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
    % make the threshold tighter to find correct inliers
    iterNum = 20000;
    thre_rsc = 2;
    sampleNum = 8;
    % RANSAC optimal result initialization
    iterOpt = 0;
    inlierNumOpt = 0;
    avrResOpt = inf;
    F_Opt = zeros(3,3);
    inliers = cell(1,iterNum);
    % run the RANSAC algorithm
    for i = 1:iterNum
        % Sampling matches (randomly choose 4)
        sampleId = randperm(num_dsc,sampleNum);
        pts1 = matches01(sampleId,:);
        pts2 = matches02(sampleId,:);
        matches_new = [pts1, pts2];

        % find an initial fitting homography based on the samples
        F = fit_fundamental(matches_new, 'normalized'); 

        % run my RANSAC function
        [inlierNum, avrRes, inliers{i}] = my_ransac(F, matches01, matches02, thre_rsc);

        if inlierNumOpt < inlierNum
            inlierNumOpt = inlierNum;
            avrResOpt = avrRes;
            iterOpt = i;
            F_Opt = F;
        end
    end
    inlierOpt = inliers{iterOpt};
    num_in = size(inlierOpt,1);
    
    %% temperory part
    imshow([I1 I2]); hold on;
    plot(inlierOpt(:,1), inlierOpt(:,2), '+r');
    plot(inlierOpt(:,3)+size(I1,2), inlierOpt(:,4), '+r');
    line([inlierOpt(:,1) inlierOpt(:,3) + size(I1,2)]', inlierOpt(:,[2 4])', 'Color', 'r');
    pause;

    %% Recalculate the optimal epipolar lines
    L = (F_Opt * [inlierOpt(:,1:2) ones(num_in,1)]')'; 
    L = L ./ repmat(sqrt(L(:,1).^2 + L(:,2).^2), 1, 3); % rescale the line
    pt_line_dist = sum(L .* [inlierOpt(:,3:4) ones(num_in,1)],2);
    closest_pt = inlierOpt(:,3:4) - L(:,1:2) .* repmat(pt_line_dist, 1, 2);

    % find endpoints of segment on epipolar line (for display purposes)
    pt1 = closest_pt - [L(:,2) -L(:,1)] * 10; % offset from the closest point is 10 pixels
    pt2 = closest_pt + [L(:,2) -L(:,1)] * 10;

    % display points and segments of corresponding epipolar lines
    clf;
    I02 = imread('library2.jpg');
    imshow(I02); hold on;
    plot(inlierOpt(:,3), inlierOpt(:,4), '+r');
    line([inlierOpt(:,3) closest_pt(:,1)]', [inlierOpt(:,4) closest_pt(:,2)]', 'Color', 'r');
    line([pt1(:,1) pt2(:,1)]', [pt1(:,2) pt2(:,2)]', 'Color', 'g');

    % use legend to report the mean squared distance btw points in both images 
    % and the corresponding epipolar lines
    legend(['inlier number = ' num2str(inlierNumOpt)  ' average Residual = ' num2str(avrResOpt)]);
    

    % report the number, the average residual, inlier locations on the plot
    figure, imagesc(I1), axis image, colormap(gray), hold on
    plot(inlierOpt(:,1),inlierOpt(:,2),'bs'), title('RANSAC inlier matches-left image'),
    legend(['inlier number = ' num2str(inlierNumOpt) ', average Residuals = ' num2str(avrResOpt)]);
    figure, imagesc(I2), axis image, colormap(gray), hold on
    plot(inlierOpt(:,3),inlierOpt(:,4),'bs'), title('RANSAC inlier matches-right image');
    legend(['inlier number = ' num2str(inlierNumOpt) ', average Residuals = ' num2str(avrResOpt)]);

    
end