%% Bonus 02:
%% Another image with threee orthogonal points and give the measurements

%% housing work
close all
clc

%% i. load images and find vainishing points

im = imread('bonus02.jpeg');

% vp_x = getVanishingPoint(im);
% pause
% vp_y = getVanishingPoint(im);
% pause
% vp_z = getVanishingPoint(im);
% pause

% results from the above functions
vp_x = [1570.4;749.3;1];
vp_y = [431.3;-7455.5;1];
vp_z = [-1218.4;747.3;1];

% plot the ground horizon line and specify its normalized parameters
figure
hold on
set(gca,'Ydir','reverse')
plot(vp_x(1), vp_x(2), 'r*');
plot(vp_y(1), vp_y(2), 'r*');
plot(vp_z(1), vp_z(2), 'r*');
axis equal
imagesc(im);

figure
imagesc(im);    
hold on;
plot([vp_x(1) vp_z(1)], [vp_x(2) vp_z(2)]);
axis image;

horiz_line = real(cross(vp_x, vp_z));
horiz_line = horiz_line/sqrt(sum(horiz_line(1:2).^2,1));


%% ii. solve for the focal length and optical center (principal point)
syms f u v;

% % first, obtain the matrix for K_inv and the middle matrix for equations
% K = [f 0 u; 0 f v; 0 0 1];
% K_inv = inv(K); 
% K_inv = [1/f, 0, -u/f; 0, 1/f, -v/f; 0, 0, 1];
% middle = transpose(K_inv) * K_inv;
% 
% middle = 
% 
% [  1/f^2,      0,                -u/f^2]
% [      0,  1/f^2,                -v/f^2]
% [ -u/f^2, -v/f^2, u^2/f^2 + v^2/f^2 + 1]

% then hard code the equations by implementing the items of matrix "middle"
[Sf, Su, Sv] = solve(...
    vp_y(1)*(vp_x(1)*1/f^2 - u/f^2) + vp_y(2)*(vp_x(2)/f^2 - v/f^2) + (vp_x(1)*(-u/f^2) + vp_x(2)*(-v/f^2)+u^2/f^2 +v^2/f^2+1) == 0,...
    vp_y(1)*(vp_z(1)*1/f^2 - u/f^2) + vp_y(2)*(vp_z(2)/f^2 - v/f^2) + (vp_z(1)*(-u/f^2) + vp_z(2)*(-v/f^2)+u^2/f^2 +v^2/f^2+1) == 0,...
    vp_z(1)*(vp_x(1)*1/f^2 - u/f^2) + vp_z(2)*(vp_x(2)/f^2 - v/f^2) + (vp_x(1)*(-u/f^2) + vp_x(2)*(-v/f^2)+u^2/f^2 +v^2/f^2+1) == 0);

% get matrix K
f = double(Sf(1));
u = double(Su(1));
v = double(Sv(1));
K = [f 0 u; 0 f v; 0 0 1];

%% iii. compute the rotation matrix for the camera
% setting the vertical vanishing point as the Y-direction, 
% the right-most vanishing point as the X-direction, 
% and the left-most vanishing point as the Z-direction.
r_x = K \ vp_x;
r_y = K \ vp_y;
r_z = K \ vp_z;

% normalize the r_x, r_y, r_z
r_x = r_x / sqrt(sum(r_x.^2,1));
r_y = r_y / sqrt(sum(r_y.^2,1));
r_z = r_z / sqrt(sum(r_z.^2,1));

R = [r_x r_y r_z];

%% iv. estimate the length of the ZJU blue part of Empire State Building
% by checking the height data in "height_data_of_EmpireStateBuilding.jpg",
% assume we know one side of the length between 72-80th floor is 108 feet
% pick up t0 and b0 by imshow command
b0 = [469 428 1];
t0 = [469 350 1];
% pick three top reference point for ZJU blue part
r_zju = [392 432 1; 392 14 1];
% assume we know one side of the length between 72-80th floor is 108 feet
H = 108; 

% the ZJU blue part of Empire State Building
line_zju1 = real(cross(b0',r_zju(1,:)'));
v_zju = real(cross(line_zju1, horiz_line));
v_zju = v_zju'/v_zju(3);

line_zju2 = real(cross(v_zju', t0'));
measure_line = real(cross(r_zju(1,:)', r_zju(2,:)'));
t_zju = real(cross(line_zju2', measure_line'));
t_zju = t_zju'/t_zju(3);

height_zju = H*sqrt(sumsqr(r_zju(1,:)-r_zju(2,:)))*sqrt(sumsqr(vp_z'-t_zju'))/(sqrt(sumsqr(t_zju'-r_zju(1,:)))*sqrt(sumsqr(vp_z'-r_zju(2,:))));

% draw pictures
figure();
imagesc(im);
hold on;
plot([vp_z(1) vp_x(1)], [vp_z(2) vp_x(2)]);
plot([v_zju(1) b0(1)], [v_zju(2) b0(2)], 'r');
plot([t0(1) b0(1)], [t0(2) b0(2)], 'r');
plot([v_zju(1) t0(1)], [v_zju(2) t0(2)], 'r');
plot([v_zju(1) t_zju(1)], [v_zju(2) t_zju(2)], 'g');
plot([r_zju(1,1) t_zju(1)], [r_zju(2,2) t_zju(2)], 'g');
plot([r_zju(1,1) v_zju(1)], [r_zju(1,2) v_zju(2)], 'g');
plot([r_zju(1,1) r_zju(2,1)], [r_zju(1,2) r_zju(2,2)], 'y');
axis equal;
axis image;

%{

results for the height of ZJU blue part of Empire State Building:
height_zju = 549.0837;

Checking the accurate answer by looking at the height data in "height_data_of_EmpireStateBuilding.jpg"
height_zju_accurate = 404+17+59+108 = 588 feet; 

They are relatively close!

%}
