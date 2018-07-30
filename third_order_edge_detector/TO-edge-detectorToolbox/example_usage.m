% Using the Third-order edge detector

% (c) LEMS, Brown University
% Amir Tamrakar (amir_tamrakar@brown.edu)
% October 2007

close all;

%% load the image

img_filename = 'problem_im1.tif';
%img_filename = '108005.jpg';
img = imread(img_filename);

%crop image
%img = img(100:200, 200:300, :);

%% call the third-order edge detector

sigma = 1;
n = 1; % interpolate 2-fold
threshold = 5;
QaD_flag = 1;

[TO_edgemap, gen_edgemap] = third_order_edge_detector(img, sigma, n, threshold, QaD_flag);

%% display the edge map (compare against the generic edge map)

figure;
%imshow(img, [0 255]);
imagesc(img, [0 255]); colormap gray; axis equal; axis off;
% display edges computed by the gradient operator
disp_edg(gen_edgemap, 'r');
% display edges whose orientations have been computed by the third-order operator
disp_edg(TO_edgemap, 'g');

