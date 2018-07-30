clear;

dst_path = 'data/';
bmap = zeros(400,300);

% %% Ellipses
% image_name = 'ellipse.jpg';
% f1 = [20 40; 176 38; 40 131; 50 220; 227 150; 33 310];
% f2 = [80 40; 263 67; 177 131; 144 200; 232 353; 152 338];
% e = [0.8; 0.6; 0.95; 0.75; 0.85; 0.55];
% figure(1)
% for i = 1:size(f1,1)
%     bmap = create_ellipse_on_binary_map( f1(i,:), f2(i,:), e(i), bmap );
%     imshow(bmap,'border','tight');
%     plot_ellipse( f1(i,:), f2(i,:), e(i) );
% end

%% Plane curve 1
image_name = 'plane_curve1.jpg';
dx = 160;
dy = 200;
amp = 50;
x = @(t) amp*(sin(2.^t) - 1.7).*cos(t)+dx;
y = @(t) amp*(sin(2.^t) - 1.7).*sin(t)+dy;
theta = linspace(0,2*pi,300);
fxy = @(x,y) (x-dx).^2 + (y-dy).^2 - (amp*(sin(2.^(atan2(y-dy,x-dx)+pi))-1.7)).^2;
figure(2)
bmap = create_shape_on_binary_map( fxy, bmap );
imshow(bmap,'border','tight'); hold on;
plot( x(theta), y(theta), 'g' ); hold off;

%% post-process
bmap = 255*uint8(bmap);
gauss_filter = fspecial('gaussian',11,2);
bmap = imfilter(bmap,gauss_filter);
imshow(bmap,'border','tight');
imwrite(bmap,[dst_path image_name])