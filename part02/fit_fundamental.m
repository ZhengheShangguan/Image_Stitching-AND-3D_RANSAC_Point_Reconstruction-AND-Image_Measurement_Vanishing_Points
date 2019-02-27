function F = fit_fundamental(matches, method)
% this function is for calculating the fundermental matrix based on
% different choices.

switch method
    case 'unnormalized'
        %% Unnormalized eight-point algorithm
        % formulate the data matrix and the right side vector
        N = size(matches,1);
        pt01 = [matches(:,1:2),ones(N,1)]';
        pt02 = [matches(:,3:4),ones(N,1)]';
        data = [repmat(pt02(1,:)',1,3).*pt01', repmat(pt02(2,:)',1,3).*pt01', pt01(1:3,:)'];

        % solve homogeneous linear system using more matches
        [eigvec,~] = eig(data'*data);
        F = reshape(eigvec(:,1),[3,3])';

        % enforce rank-2 constraint by taking SVD of F and throwing out the
        % smallest singualr value
        [U,S,V] = svd(F);
        S(3,3) = 0;
        F = U*S*(V');
        
    case 'normalized'
        %% Normalized eight-point algorithm
        % get normalization matrix first
        % center the image data at the origin, and scale it to keep the
        % mean squared distance btw the origin and the data points is 2
        % pixels
        N = size(matches,1);
        pt01 = [matches(:,1:2),ones(N,1)]';
        pt02 = [matches(:,3:4),ones(N,1)]';
        center01 = mean(pt01,2);
        center02 = mean(pt02,2);
        dists01 = sqrt(sum((pt01 - repmat(center01,1,size(pt01,2))).^2,1));
        dists02 = sqrt(sum((pt02 - repmat(center02,1,size(pt02,2))).^2,1));
        mean_dists01 = mean(dists01);
        mean_dists02 = mean(dists02);
        norm_mat01 = [2/mean_dists01 0 -2/mean_dists01*center01(1);...
                      0 2/mean_dists01 -2/mean_dists01*center01(2);...
                      0 0 1];
        norm_mat02 = [2/mean_dists02 0 -2/mean_dists02*center02(1);...
                      0 2/mean_dists02 -2/mean_dists02*center02(2);...
                      0 0 1];
                
        % formulate the data matrix by using normalized matches
        pt01n = norm_mat01 * pt01;
        pt02n = norm_mat02 * pt02;
        data = [repmat(pt02n(1,:)',1,3).*pt01n', repmat(pt02n(2,:)',1,3).*pt01n', pt01n(1:3,:)'];

        % solve homogeneous linear system using more matches
        [eigvec,~] = eig(data'*data);
        F = reshape(eigvec(:,1),[3,3])';

        % enforce rank-2 constraint by taking SVD of F and throwing out the
        % smallest singualr value
        [U,S,V] = svd(F);
        S(3,3) = 0;
        F = U*S*V';
        
        % transform the fundamental matrix back to the original units
        F= norm_mat02' * F * norm_mat01;
        
        
    otherwise
        disp('Unknown method. Please recheck your insert.')

end

end