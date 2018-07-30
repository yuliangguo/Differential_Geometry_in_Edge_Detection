close all;
clear;
clc;

addpath( genpath('third_order_demo/'))
addpath( genpath('toolbox/'))


src_path = 'data/';
dst_path = src_path;
files = dir([src_path '*.jpg']);
thresh = 2;
sigma = 2;
% n = 1;
n = 1; % no subpixel

opts.nrad = 3.5;
opts.gap = 1.5;
opts.dx = 0.4;
opts.dt = 15;
opts.token_len = 1;
opts.max_k = 0.3;
opts.cvlet_style = 3;
opts.max_size_to_group = 7;
opts.output_type = 0;


for i = 1:length(files)
    img = imread([src_path files(i).name]);
    gauss_filter = fspecial('gaussian',11,2);
    img = imfilter(img,gauss_filter);
    imshow(img,'border','tight');
    
    opts.w = size(img,2);
    opts.h = size(img,1);
    [TO_edge_map, gen_edge_map] = third_order_edge_detector(img, sigma, n, thresh, 1);
    imshow(img,'border','tight');
    disp_edg(TO_edge_map,'g');
%     [TO_edge_map, ~] = multi_spect_TO_edge_detector(img, n, sigma, thresh, 'Lab');
%     TO_edge_map = load_edg([dst_path files(i).name(1:end-4) '.edg']);
    % TO_edge_map = [[0.5:8.5]', [0.5:8.5]', -ones(9,1)*pi/4, ones(9,1)*10];
    % TO_edge_map = [TO_edge_map; 5.1, 5.1, -pi/4, 10; 6.8, 6.8, -pi/4, 10];
    
    % NOTE
    [chain, info] = form_curvelet_mex(TO_edge_map, size(img,1), size(img,2),...
     opts.nrad, opts.gap, opts.dx, opts.dt/180*pi, opts.token_len, opts.max_k,...
     opts.cvlet_style, opts.max_size_to_group, opts.output_type);
    % !!!! the 2nd & 3rd colunms of info are kMAX & kMIN
    
    save_edg([dst_path files(i).name(1:end-4) '.edg'], TO_edge_map,[size(img,2) size(img,1)])
    save_cvlet([dst_path files(i).name(1:end-4) '.cvlet'], chain, info, TO_edge_map, opts);
%     krange = info(:,2:3);
%     save([dst_path files(i).name(1:end-4) '_krange.mat'], 'krange');
    
%     % update the edges
%     TO_edge_map2 = update_edges_via_curvelets(TO_edge_map, chain,info);
%     save_edg([dst_path files(i).name(1:end-4) '_new.edg'], TO_edge_map2,[size(img,2) size(img,1)])
% %     
%     [ cvfrag ] = convert_cfrag( chain, TO_edge_map2 );
% %     figure(2);
%     imshow(img,'border','tight');
%     hold on;
%     draw_contours(cvfrag(randperm(length(cvfrag),100)));
%     hold off;
% %     save_edg(filename, edg, dim)
end