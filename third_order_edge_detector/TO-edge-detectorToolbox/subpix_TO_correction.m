function [ TO_gradmap, TO_orientation, edginfo ] = subpix_TO_correction( gradmap,theta,thresh,mask,sigma,n)
% Use NMS to detect the edgels - this NMS routine ouputs subpixel edgel tokens
if nargin<4
    mask = gradmap>thresh;
else
    mask = mask.*(gradmap>thresh);
end

if nargin<5
    sigma = 0.7;
end
if nargin<6
    n = 0;
end
margin = 4;
[r, c] = size(gradmap);
if n
    gradmap_i = imresize(gradmap,[r*2^n+1, c*2^n+1]);
    theta_i = imresize(theta,[r*2^n+1, c*2^n+1]);
    % edgemap_i = imresize(edgemap,[r*2^n+1, c*2^n+1],'nearest');
    gradmap_i = gradmap_i(1:(r-1)*2^n+1, 1:(c-1)*2^n+1);
    theta_i = theta_i(1:(r-1)*2^n+1, 1:(c-1)*2^n+1);
    % edgemap_i = edgemap_i(1:(r-1)*2^n+1, 1:(c-1)*2^n+1);
else
    gradmap_i = gradmap;
    theta_i = theta;
end

dirx = gradmap.*cos(mod(theta+pi/2, pi));
diry = gradmap.*sin(mod(theta+pi/2, pi));
[subpix_x, subpix_y, subpix_dir_x, subpix_dir_y, subpix_grad_x, subpix_grad_y]...
    = NMS_token(dirx, diry,gradmap, mask, margin);

subpix_x = subpix_x*(2^n);
subpix_y = subpix_y*(2^n);

% magnitude of the gradient at the maxima
mag_e = sqrt(subpix_grad_x.^2+subpix_grad_y.^2);

Hx = gradmap_i.*cos(mod(theta_i+pi/2, pi));
Hy = gradmap_i.*sin(mod(theta_i+pi/2, pi));

Px   = resampled_filter_2d(gradmap, @Gx_2d_op, sigma, n);
Py   = resampled_filter_2d(gradmap, @Gy_2d_op, sigma, n);
Pxx  = resampled_filter_2d(gradmap, @Gxx_2d_op, sigma, n);
Pxy  = resampled_filter_2d(gradmap, @Gxy_2d_op, sigma, n);
Pyy  = resampled_filter_2d(gradmap, @Gyy_2d_op, sigma, n);

Hxx  = resampled_filter_2d(dirx, @Gx_2d_op, sigma, n);
Hyy  = resampled_filter_2d(diry, @Gy_2d_op, sigma, n);
Hxy  = resampled_filter_2d(diry, @Gx_2d_op, sigma, n);

[yd, xd] = size(Hx);
[xx,yy] = meshgrid(1:xd,1:yd);
% grad_e = interp2(xx, yy, gradmap_i, subpix_x, subpix_y, 'cubic');
Px_e   = interp2(xx, yy, Px,  subpix_x, subpix_y, 'cubic');
Py_e   = interp2(xx, yy, Py,  subpix_x, subpix_y, 'cubic');
Pxx_e  = interp2(xx, yy, Pxx,  subpix_x, subpix_y, 'cubic');
Pyy_e  = interp2(xx, yy, Pyy,  subpix_x, subpix_y, 'cubic');
Pxy_e  = interp2(xx, yy, Pxy,  subpix_x, subpix_y, 'cubic');
Hx_e  = interp2(xx, yy, Hx, subpix_x, subpix_y, 'cubic');
Hy_e  = interp2(xx, yy, Hy, subpix_x, subpix_y, 'cubic');
Hxx_e = interp2(xx, yy, Hxx, subpix_x, subpix_y, 'cubic');
Hxy_e = interp2(xx, yy, Hxy, subpix_x, subpix_y, 'cubic');
Hyy_e = interp2(xx, yy, Hyy, subpix_x, subpix_y, 'cubic');

% compute [Fx, Fy] at zero crossings (second order equations on P = |hist-grad I|)
Fx_e = Px_e.*Hxx_e + Pxx_e.*Hx_e + Py_e.*Hxy_e + Pxy_e.*Hy_e;
Fy_e = Px_e.*Hxy_e + Pxy_e.*Hx_e + Py_e.*Hyy_e + Pyy_e.*Hy_e;
% Fx_e = Px_e.*Hxx_e + Pxx_e.*Hx_e + Py_e.*Hxy_e + Pxy_e.*Hy_e - Px_e.*(Px_e.*Hx_e + Py_e.*Hy_e)./grad_e;
% Fy_e = Px_e.*Hxy_e + Pxy_e.*Hx_e + Py_e.*Hyy_e + Pyy_e.*Hy_e - Py_e.*(Px_e.*Hx_e + Py_e.*Hy_e)./grad_e;

%normalize
F_mag = sqrt(Fx_e.^2 + Fy_e.^2);
Fx_e = Fx_e./F_mag;
Fy_e = Fy_e./F_mag;

% Edge Direction (tangent to the level set) is orthogonal to the gradient
subpix_dir_x2 = (-Fy_e);
subpix_dir_y2 = (Fx_e);

% construct the edge map
edginfo = [subpix_x'/(2^n)-1 subpix_y'/(2^n)-1 atan2(subpix_dir_y2', subpix_dir_x2') mag_e']; 
TO_gradmap = zeros(r,c);
TO_orientation = zeros(r,c);
for q = 1:size(edginfo, 1)
    x = min([max([1,round(edginfo(q, 1)+1)]),c]);
    y = min([max([1,round(edginfo(q, 2)+1)]),r]);
    TO_gradmap(y,x) = edginfo(q, 4);
    TO_orientation(y,x) = edginfo(q, 3);
end
        
% debug only
% figure;
% disp_edg(edginfo, 'b');
% axis([100 180 100 180]);

end

