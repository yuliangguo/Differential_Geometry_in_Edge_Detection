function binary_map = create_ellipse_on_binary_map( f1, f2, e, background )
% Input:
%       f1, f2: [1x2] focal points;
%       e: the eccentricity of an ellipse, f/a;
%       background: binary image background
[h,w] = size(background);
[xx,yy] = meshgrid(1:w,1:h);
a = 1/2*sqrt((f2(1)-f1(1))^2+(f2(2)-f1(2))^2);
b = a*sqrt(1-e^2);
center = (f1+f2)/2;
theta = -atan2(f2(2)-f1(2),f2(1)-f1(1));
rotation = [cos(theta) -sin(theta); sin(theta) cos(theta)];
trans_pos = rotation*[xx(:)'-center(1);yy(:)'-center(2)];
ind = (trans_pos(1,:)/a).^2+(trans_pos(2,:)/b).^2-1<0;
background(ind) = 1;
binary_map = background;
end