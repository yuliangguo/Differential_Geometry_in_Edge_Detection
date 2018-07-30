addpath(genpath('../../toolbox/'))
src_path = 'data/';
file_name = 'plane_curve1.cem';%'circle_fill_blur.cem';
img_name = 'plane_curve1.jpg';
img = imread([src_path img_name]);
cem = load_contours([src_path file_name]);
figure(1);
imshow(img,'border','tight');
% disp_edg(cat(1,cem{2}{:}),'g');
% show the curves with curve id to find out the curvature GT;
hold on;
colormap = hsv(length(cem{2}));
for i = 1:length(cem{2})
    plot(cem{2}{i}(:,1)+1,cem{2}{i}(:,2)+1,'color',colormap(i,:));
    mid_id = ceil(size(cem{2}{i}, 1)/2);
    text([cem{2}{i}(mid_id,1)], ...
        [cem{2}{i}(mid_id,2)], ...
        num2str(i),'color',colormap(i,:), 'FontSize',14);
end
hold off;

