function H = findHomo(pts1, pts2, sampleNum)

%% Calculating a fitting homography
% 1. Define matrix A
A = zeros(sampleNum*2,9);
A(1:2:2*sampleNum,4:5) = pts1;
A(1:2:2*sampleNum,6) = 1;
A(2:2:2*sampleNum,1:2) = pts1;
A(2:2:2*sampleNum,3) = 1;
x1 = pts1(:,1);
y1 = pts1(:,2);
x2 = pts2(:,1);
y2 = pts2(:,2);
A(1:2:2*sampleNum,7) = -y2.*x1;
A(1:2:2*sampleNum,8) = -y2.*y1;
A(1:2:2*sampleNum,9) = -y2;
A(2:2:2*sampleNum,7) = -x2.*x1;
A(2:2:2*sampleNum,8) = -x2.*y1;
A(2:2:2*sampleNum,9) = -x2;
% 2. Calculate H
[eigvec,~] = eig(A'*A);
H = reshape(eigvec(:,1),[3,3])';
H = H/H(end); % make H(3,3)=1 thus cancel lambda factor in the formula

end