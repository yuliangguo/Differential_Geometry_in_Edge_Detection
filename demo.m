clear all; close all;
addpath(genpath('./'));

input_name = '2018';
dst_path = './';
img = imread([input_name, '.jpg']);
opts.w = size(img,2);
opts.h = size(img,1);

%% 1. Third-Order Edge Detection Gray
disp('third-order edge detection: gray')
% edge detector parameters (variables)
thresh = 2;
sigma = 2;
n = 1;
tic;
[TO_edges, gen_edges] = third_order_edge_detector(img, sigma, n, thresh, 1);
toc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% outputs:
% TO_edges = [Subpixel_X Subpixel_Y Orientation Confidence]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% 2. %%%%%%%%%%%%%%%%%  Curvature estimation
disp('curvature estimation')
% n = 0; % no subpixel
% Curvature Estimation parameters (Fixed)
opts.nrad = 3.5;
opts.gap = 1.5;
opts.dx = 0.4;
opts.dt = 15;
opts.token_len = 1;
opts.max_k = 0.3;
opts.cvlet_style = 3;
opts.max_size_to_group = 7;
opts.output_type = 0;
% perform Curvature Estimation
[chain, info] = form_curvelet_mex(TO_edges, size(img,1), size(img,2),...
opts.nrad, opts.gap, opts.dx, opts.dt/180*pi, opts.token_len, opts.max_k,...
opts.cvlet_style, opts.max_size_to_group, opts.output_type);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% outputs:
% chain: edge id chains of curves
% info: [isForward ref_pt.x() ref_pt.y() ref_theta pt.x() pt.y() theta k length property]
% use pt.x(), pt.y(), k to output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % in case you want to output edges and curvals in files
% save_edg([dst_path input_name '.edg'], TO_edges,[size(img,2) size(img,1)])
% save_cvlet([dst_path input_name '.cvlet'], chain, info, TO_edges, opts);

%% 3. visualization
disp('display curvature map')
k_map = zeros(opts.h, opts.w);
info(:,5) = min(info(:,5), opts.w-1);
info(:,6) = min(info(:,6), opts.h-1);

for i = 1:size(info,1)
    k_map(round(info(i,6))+1, round(info(i,5))+1) = max(k_map(round(info(i,6))+1, round(info(i,5))+1), abs(info(i,8)));
end
figure(1);
imagesc(k_map);
axis image

disp('display edges')
figure(2)
imshow(img, 'border', 'tight')
% This visulization can be slow when there are many edges to show. Users are recommend to zoom
% in to look the subpixel location and accurate orientation of detected edges
disp_edg(TO_edges);

%% 4. Third-Order Edge Detection Multi-channel
disp('thrid-order edge detection: color')
tic;
channel_opt = {'rgb','Lab','HSV','oppo'};
TO_edges = multi_spect_TO_edge_detector(img, n, sigma, thresh, channel_opt{2}, 0);
toc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% outputs:
% TO_edges2 = [Subpixel_X Subpixel_Y Orientation Confidence]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 5. %%%%%%%%%%%%%%%%% Third-Party Boundary Detector 
disp('Structured edge detector with third-order edge orientation and subpixel edge localization')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this part is based on
% "Fast Edge Detection Using Structured Forests"
% Piotr Dollar and C. Lawrence Zitnick
% IEEE Transactions on Pattern Analysis and Machine Intelligence 2015
% The user need set up their code package and run the following code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
th = 0.1; % threshold for boundary map to localize edges
sigma = 1; % NMS kernal size
% Original SE edge detection 
[E,O,inds,segs] = edgesDetect( I, model); 
% Third-order orientation estimation and subpixel localization
[ TO_edgemap, TO_orientation, TO_edges] = subpix_TO_correction( E,O,th,sigma);
