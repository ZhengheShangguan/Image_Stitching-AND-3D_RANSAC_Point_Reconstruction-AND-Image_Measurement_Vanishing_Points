function [inlierNum, avrRes, inliers] = my_ransac(F, matches01, matches02, thre_rsc)

%% Find out inliers with thre_rsc


% transform points from the first image to get epipolar lines in the second image
num_dsc = size(matches01, 1);
L = (F * [matches01 ones(num_dsc,1)]')'; 

% find points on epipolar lines L closest to matches(:,3:4)
L = L ./ repmat(sqrt(L(:,1).^2 + L(:,2).^2), 1, 3); % rescale the line
pt_line_dist = sum(L .* [matches02 ones(num_dsc,1)],2);
closest_pt = matches02 - L(:,1:2) .* repmat(pt_line_dist, 1, 2);

% the mean squared distance btw points in both images and the corresponding epipolar lines
dist = sqrt(sum((closest_pt - matches02).^2,2));
% average Residual
avrRes = mean(dist);


% find the inliers according to their distances, calculate the whole residual
index = find(dist<thre_rsc);
inliers = zeros(length(index),4);
for i = 1:length(index)
    inliers(i,:) = [matches01(index(i),1:2), matches02(index(i),1:2)];
end
% calculate inlierNum and avrRes
inlierNum = size(inliers,1);

end

