clc;
clear;
addpath(genpath('./'))
path='./';
savepath = './';
saveedg = {'rgb','Lab','HSV','oppo'};
saveimg = {'rgb_img','Lab_img','HSV_img','oppo_img'};
% ext='.jpg';
% imgpath = LoadImages(path,ext);
interp=1;
sigma=2;
TH=1;
for j=1:4
    dirpath = [savepath '/' saveedg{j} '_T1/'];
    if ~isdir(dirpath)
        mkdir(dirpath);
    end
    dirpathimg = [savepath '/' saveimg{j} '_T1/'];
    if ~isdir(dirpathimg)
        mkdir(dirpathimg);
    end
        img=imread('2018.jpg');
        [y,x,z]=size(img);
        dim = [x,y];
        TO_edge=multi_spect_TO_edge_detector(img, interp, sigma, TH, saveedg{j},0);
    %     [edgex, edgey, orientation]=gray_TO_edge_detector(img, 'r');
        TO_edge(:,1:2) = TO_edge(:,1:2) - 1;
        name = '2018';
%         [~,name,~] = fileparts(imgpath{i});
        saveedgname = [dirpath 'edge_' name '.edg'];
        save_edg(saveedgname, TO_edge, dim);
%         pause(1);
        saveimgname = [dirpathimg 'edge_' name];
        saveas(gcf, saveimgname, 'jpg');
        close all;
end