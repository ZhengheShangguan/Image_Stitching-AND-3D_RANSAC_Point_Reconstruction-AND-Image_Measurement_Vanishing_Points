function [inlierNum, avrRes, inliers] = my_ransac(H, pts1, pts2, thre_rsc)

%% Find out inliers with thre_rsc
% reprojection and calculating the distances
n = size(pts1,1);
pts3 = H*[pts1,ones(n,1)]';
pts3 = pts3(1:2,:)./repmat(pts3(3,:),2,1);
dist = sum((pts2'-pts3).^2,1);
% find the inliers according to their distances, calculate the whole residual
index = find(dist<thre_rsc);
inliers = zeros(length(index),4);
avrRes = 0;
for i = 1:length(index)
    inliers(i,:) = [pts1(index(i),1:2), pts2(index(i),1:2)];
    avrRes = avrRes + dist(index(i));
end
% calculate inlierNum and avrRes
inlierNum = size(inliers,1);
avrRes = avrRes/inlierNum;

end

