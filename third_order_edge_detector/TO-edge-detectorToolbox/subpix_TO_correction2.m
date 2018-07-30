function [ TO_gradmap, TO_orientation, edginfo ] = subpix_TO_correction2( gradmap,theta,thresh,sigma,n)
% Use NMS to detect the edgels - this NMS routine ouputs subpixel edgel tokens
if nargin<4
    sigma = 0.7;
end
if nargin<5
    n = 0;
end
margin = 4;
if n
    [r, c] = size(gradmap);
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

dirx = sin(theta);
diry = -cos(theta);
[subpix_x, subpix_y, subpix_dir_x, subpix_dir_y, subpix_grad_x, subpix_grad_y]...
    = NMS_token(dirx, diry,gradmap, gradmap>thresh, margin);

subpix_x = subpix_x*(2^n);
subpix_y = subpix_y*(2^n);

% magnitude of the gradient at the maxima
mag_e = sqrt(subpix_grad_x.^2+subpix_grad_y.^2);

sintheta = sin(theta_i);
costheta = cos(theta_i);

Px   = resampled_filter_2d(gradmap, @Gx_2d_op, sigma, n);
Py   = resampled_filter_2d(gradmap, @Gy_2d_op, sigma, n);
Pxx  = resampled_filter_2d(gradmap, @Gxx_2d_op, sigma, n);
Pxy  = resampled_filter_2d(gradmap, @Gxy_2d_op, sigma, n);
Pyy  = resampled_filter_2d(gradmap, @Gyy_2d_op, sigma, n);

thetax  = resampled_filter_2d(theta, @Gx_2d_op, sigma, n);
thetay  = resampled_filter_2d(theta, @Gy_2d_op, sigma, n);

[yd, xd] = size(sintheta);
[xx,yy] = meshgrid(1:xd,1:yd);
% grad_e = interp2(xx, yy, gradmap_i, subpix_x, subpix_y, 'cubic');
Px_e   = interp2(xx, yy, Px,  subpix_x, subpix_y, 'cubic');
Py_e   = interp2(xx, yy, Py,  subpix_x, subpix_y, 'cubic');
Pxx_e  = interp2(xx, yy, Pxx,  subpix_x, subpix_y, 'cubic');
Pyy_e  = interp2(xx, yy, Pyy,  subpix_x, subpix_y, 'cubic');
Pxy_e  = interp2(xx, yy, Pxy,  subpix_x, subpix_y, 'cubic');
sintheta_e  = interp2(xx, yy, sintheta, subpix_x, subpix_y, 'cubic');
costheta_e  = interp2(xx, yy, costheta, subpix_x, subpix_y, 'cubic');
thetax_e = interp2(xx, yy, thetax, subpix_x, subpix_y, 'cubic');
thetay_e = interp2(xx, yy, thetay, subpix_x, subpix_y, 'cubic');

% compute [Fx, Fy] at zero crossings (second order equations on P = |hist-grad I|)
Fx_e = Px_e.*costheta_e.*thetax_e + Pxx_e.*sintheta_e + Py_e.*sintheta_e.*thetax_e...
    - Pxy_e.*costheta_e;
Fy_e = Px_e.*costheta_e.*thetay_e + Pxy_e.*sintheta_e + Py_e.*sintheta_e.*thetay_e...
    - Pyy_e.*costheta_e;

%normalize
F_mag = sqrt(Fx_e.^2 + Fy_e.^2);
Fx_e = Fx_e./F_mag;
Fy_e = Fy_e./F_mag;

% Edge Direction (tangent to the level set) is orthogonal to the gradient
subpix_dir_x2 = (-Fy_e);
subpix_dir_y2 = (Fx_e);

% construct the edge map
edginfo = [subpix_x'/(2^n)-1 subpix_y'/(2^n)-1 atan2(subpix_dir_y2', subpix_dir_x2') mag_e']; 
TO_gradmap = zeros(size(gradmap));
TO_orientation = zeros(size(gradmap));
for q = 1:size(edginfo, 1)
    x = round(edginfo(q, 1)+1);
    y = round(edginfo(q, 2)+1);
    TO_gradmap(y,x) = edginfo(q, 4);
    TO_orientation(y,x) = edginfo(q, 3);
end
        
% debug only
% figure;
% disp_edg(edginfo, 'b');
% axis([100 180 100 180]);

end

