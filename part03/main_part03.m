%% 
%% housing work
close all
clc

%% i. load images and find vainishing points

im = imread('CSL.jpg');

% vp_x = getVanishingPoint(im);
% pause
% vp_y = getVanishingPoint(im);
% pause
% vp_z = getVanishingPoint(im);
% pause

% results from the above functions
vp_x = [1344.1;228.8;1];
vp_y = [-239.1;213.7;1];
vp_z = [519.4;6742.5;1];

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
plot([vp_x(1) vp_y(1)], [vp_x(2) vp_y(2)]);
axis image;

horiz_line = real(cross(vp_x, vp_y));
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
tmp = vp_z;
vp_z = vp_y;
vp_y = tmp;
r_x = K \ vp_x;
r_y = K \ vp_y;
r_z = K \ vp_z;

% normalize the r_x, r_y, r_z
r_x = r_x / sqrt(sum(r_x.^2,1));
r_y = r_y / sqrt(sum(r_y.^2,1));
r_z = r_z / sqrt(sum(r_z.^2,1));

R = [r_x r_y r_z];

%% iv. estimate the heights of
% pick up t0 and b0 by imshow command
b0 = [628 510 1];
t0 = [628 467 1];
% pick three top reference point for CSL, Spike and lamp post
r_csl = [602 307 1; 602 53 1];
r_spike = [603 467 1; 603 190 1];
r_lamp = [294 522 1; 294 390 1];
% assuming the person is 5ft 6in (66 inch)/ 6ft (72 inch) tall
H = 66; % 72;

% (a) the CSL building
line_csl1 = real(cross(b0',r_csl(1,:)'));
v_csl = real(cross(line_csl1, horiz_line));
v_csl = v_csl'/v_csl(3);

line_csl2 = real(cross(v_csl', t0'));
measure_line = real(cross(r_csl(1,:)', r_csl(2,:)'));
t_csl = real(cross(line_csl2', measure_line'));
t_csl = t_csl'/t_csl(3);

height_csl = H*sqrt(sumsqr(r_csl(1,:)-r_csl(2,:)))*sqrt(sumsqr(vp_z'-t_csl'))/(sqrt(sumsqr(t_csl'-r_csl(1,:)))*sqrt(sumsqr(vp_z'-r_csl(2,:))));

% draw pictures
figure();
imagesc(im);
hold on;
plot([vp_z(1) vp_x(1)], [vp_z(2) vp_x(2)]);
plot([v_csl(1) b0(1)], [v_csl(2) b0(2)], 'r');
plot([t0(1) b0(1)], [t0(2) b0(2)], 'r');
plot([v_csl(1) t0(1)], [v_csl(2) t0(2)], 'r');
plot([v_csl(1) t_csl(1)], [v_csl(2) t_csl(2)], 'g');
plot([r_csl(1,1) t_csl(1)], [r_csl(2,2) t_csl(2)], 'g');
plot([r_csl(1,1) v_csl(1)], [r_csl(1,2) v_csl(2)], 'g');
plot([r_csl(1,1) r_csl(2,1)], [r_csl(1,2) r_csl(2,2)], 'y');
axis equal;
axis image;

% (b) the spike statue
line_spike1 = real(cross(b0',r_spike(1,:)'));
v_spike = real(cross(line_spike1, horiz_line));
v_spike = v_spike'/v_spike(3);

line_spike2 = real(cross(v_spike', t0'));
measure_line = real(cross(r_spike(1,:)', r_spike(2,:)'));
t_spike = real(cross(line_spike2', measure_line'));
t_spike = t_spike'/t_spike(3);

height_spike = H*sqrt(sumsqr(r_spike(1,:)-r_spike(2,:)))*sqrt(sumsqr(vp_z'-t_spike'))/(sqrt(sumsqr(t_spike'-r_spike(1,:)))*sqrt(sumsqr(vp_z'-r_spike(2,:))));

% draw pictures
figure();
imagesc(im);
hold on;
plot([vp_z(1) vp_x(1)], [vp_z(2) vp_x(2)]);
plot([v_spike(1) b0(1)], [v_spike(2) b0(2)], 'r');
plot([t0(1) b0(1)], [t0(2) b0(2)], 'r');
plot([v_spike(1) t0(1)], [v_spike(2) t0(2)], 'r');
plot([v_spike(1) t_spike(1)], [v_spike(2) t_spike(2)], 'g');
plot([r_spike(1,1) t_spike(1)], [r_spike(2,2) t_spike(2)], 'g');
plot([r_spike(1,1) v_spike(1)], [r_spike(1,2) v_spike(2)], 'g');
plot([r_spike(1,1) r_spike(2,1)], [r_spike(1,2) r_spike(2,2)], 'y');
axis equal;
axis image;

% (c) the lamp posts
line_lamp1 = real(cross(b0',r_lamp(1,:)'));
v_lamp = real(cross(line_lamp1, horiz_line));
v_lamp = v_lamp'/v_lamp(3);

line_lamp2 = real(cross(v_lamp', t0'));
measure_line = real(cross(r_lamp(1,:)', r_lamp(2,:)'));
t_lamp = real(cross(line_lamp2', measure_line'));
t_lamp = t_lamp'/t_lamp(3);

height_lamp = H*sqrt(sumsqr(r_lamp(1,:)-r_lamp(2,:)))*sqrt(sumsqr(vp_z'-t_lamp'))/(sqrt(sumsqr(t_lamp'-r_lamp(1,:)))*sqrt(sumsqr(vp_z'-r_lamp(2,:))));

% draw pictures
figure();
imagesc(im);
hold on;
plot([vp_z(1) vp_x(1)], [vp_z(2) vp_x(2)]);
plot([v_lamp(1) b0(1)], [v_lamp(2) b0(2)], 'r');
plot([t0(1) b0(1)], [t0(2) b0(2)], 'r');
plot([v_lamp(1) t0(1)], [v_lamp(2) t0(2)], 'r');
plot([v_lamp(1) t_lamp(1)], [v_lamp(2) t_lamp(2)], 'g');
plot([r_lamp(1,1) t_lamp(1)], [r_lamp(2,2) t_lamp(2)], 'g');
plot([r_lamp(1,1) v_lamp(1)], [r_lamp(1,2) v_lamp(2)], 'g');
plot([r_lamp(1,1) r_lamp(2,1)], [r_lamp(1,2) r_lamp(2,2)], 'y');
axis equal;
axis image;

%{

results for the heights:
height_csl = 1299.3
height_spike = 515.35;
height_lamp = 203.76;

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Bonus points

%% 1. perform additional measurements
%%%%%%%%%%%%%%%%%%
% the tallest man
% the height of the top point of the front large window on CSL Building
% the height of the tree on the left side of the front door of CSL

% different objects to choose for the measurement
r_firstman = [b0; t0];
r_orange = [473 577 1; 473 520 1];
r_black = [455 565 1; 455 516 1];
r_white = [322 402 1; 322 374 1];
r_window = [487 310 1; 487 158 1];
r_tree = [364 308 1; 364 234 1];

% choose different objects from the above lists: man, window or tree
r_man = r_tree;

line_man1 = real(cross(b0',r_man(1,:)'));
v_man = real(cross(line_man1, horiz_line));
v_man = v_man'/v_man(3);

line_man2 = real(cross(v_man', t0'));
measure_line = real(cross(r_man(1,:)', r_man(2,:)'));
t_man = real(cross(line_man2', measure_line'));
t_man = t_man'/t_man(3);

height_man = H*sqrt(sumsqr(r_man(1,:)-r_man(2,:)))*sqrt(sumsqr(vp_z'-t_man'))/(sqrt(sumsqr(t_man'-r_man(1,:)))*sqrt(sumsqr(vp_z'-r_man(2,:))));

% draw pictures
figure();
imagesc(im);
hold on;
plot([vp_z(1) vp_x(1)], [vp_z(2) vp_x(2)]);
plot([v_man(1) b0(1)], [v_man(2) b0(2)], 'r');
plot([t0(1) b0(1)], [t0(2) b0(2)], 'r');
plot([v_man(1) t0(1)], [v_man(2) t0(2)], 'r');
plot([v_man(1) t_man(1)], [v_man(2) t_man(2)], 'g');
plot([r_man(1,1) t_man(1)], [r_man(2,2) t_man(2)], 'g');
plot([r_man(1,1) v_man(1)], [r_man(1,2) v_man(2)], 'g');
plot([r_man(1,1) r_man(2,1)], [r_man(1,2) r_man(2,2)], 'y');
axis equal;
axis image;

%{

1. results for the tallest man:
height_orange = 70.8;
height_black = 62.7;
height_white = 67.6;
height_firstman = 66;
Thus, the man in orange is the tallest man in the picture.

2. results for the top point of the CSL window:
height_window = 754.54;

3. results for the height of the tree on the left side of CSL front door
height_tree = 372.57;

%}

%% 2. Another image with threee orthogonal points and give the measurements

% Please see the codes in the script "bonus02.m", thanks!


