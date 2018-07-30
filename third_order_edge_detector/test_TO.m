clear;
addpath TO-edge-detectorToolbox/
addpath(genpath('../toolbox'))


% for d = 1:7
% img_path = ['/home/xiaoyan/Documents/MATLAB/Third_order/eval/camera_noise/Data/cropSET'...
%     num2str(d) '/'];
% dst_path = ['/home/xiaoyan/Documents/MATLAB/Third_order/eval/camera_noise/Data/edg_TO/SET'...
%     num2str(d) '/'];
% img_path = '/home/xiaoyan/Documents/MATLAB/Third_order/Data/VOC2007_data/VOC2007_img_scales/s1/';
% dst_path = '/home/xiaoyan/Documents/MATLAB/Third_order/Data/VOC2007_RES/VOC2007_edg_scales/TO/';
% img_path = ['/home/xiaoyan/Documents/MATLAB/Third_order/eval/gloabal_measure_of_coherence/Data/benchmark'...
%     num2str(d) '/images/'];
% dst_path = ['/home/xiaoyan/Documents/MATLAB/Third_order/eval/gloabal_measure_of_coherence/Data/benchmark'...
%     num2str(d) '/results/'];
img_path = '/home/xiaoyan/Documents/MATLAB/Third_order/Data/VOC2007_data/VOC2007_img_scales/s1/';
dst_path = '/home/xiaoyan/Documents/MATLAB/Third_order/Data/VOC2007_RES/VOC2007_edg_scales/TO/';

mkdir(dst_path);
% mkdir([dst_path 'edg_grad/']);
% mkdir([dst_path 'edg_TO/']);

% set edge detection parameter
sigma = 1.5;
n = 1; % interpolate 2-fold
threshold = 0;
QaD_flag = 1;

files=dir([img_path '*jpg']);

for j=1:length(files)
    name=files(j).name;
    imgFile=[img_path, name];
    fprintf(['file: ', num2str(j), '/', num2str(length(files)), '\n']);
    
    outedg=[dst_path, name(1:end-4) ];
%     outedg=[dst_path, 'edg_TO/', name(1:end-4) ];
%     if exist([outedg '.edg'],'file'), continue; end
    
    img = imread (imgFile);
%     [TO_edgemap, gen_edgemap] = third_order_edge_detector(img, sigma, n, threshold, QaD_flag);
    [TO_edgemap, ~] = multi_spect_TO_edge_detector(img, n, sigma, threshold, 'Lab');
%     [~, gen_edgemap] = multi_spect_TO_edge_detector(img, 0, sigma, threshold, 'Lab');
    siz = size(img);
    
    % convert to GMC coordinate
%     TO_edgemap(:,2) = siz(1)-1-TO_edgemap(:,2);
%     TO_edgemap(:,3) = -TO_edgemap(:,3);
%     save_edg([dst_path 'edg_TO/' name(1:end-4) '.edg'], TO_edgemap, [siz(2),siz(1)])
%     gen_edgemap(:,2) = siz(1)-1-gen_edgemap(:,2);
%     gen_edgemap(:,3) = -gen_edgemap(:,3);
%     save_edg([dst_path 'edg_grad/' name(1:end-4) '.edg'], gen_edgemap, [siz(2),siz(1)])
    
    save_edg([outedg '.edg'], TO_edgemap, [siz(2),siz(1)])
%     save_edg([outedg '.edg'], gen_edgemap, [siz(2),siz(1)])
end
% end