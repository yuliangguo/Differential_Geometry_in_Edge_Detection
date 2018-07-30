function [edg, tg, theta] = tg_subpixel_TO(img, thresh)
%
% compute subpixel edges from pb response and correct the orientation 
% using the third-order method
%
% opt: 0=pBG, 1:pTG, 2:pBGTG, 3:pCG, 4:pCGTG
%
% parameters
radius=[0.01 0.02 0.02 0.02];
norient=8;

    im = im2double(img);
    [h w ~] = size(im);
    idiag = norm([h w]);
    
    % compute texture gradient
    no = 6;
    ss = 1;
    ns = 2;
    sc = sqrt(2);
    el = 2;
    k = 64;
    fname = sprintf( ...
        'unitex_%.2g_%.2g_%.2g_%.2g_%.2g_%d.mat',no,ss,ns,sc,el,k);
    textonData = load(fname); % defines fb,tex,tsim
    if size(im,3)==3
        im=rgb2gray(im);
    end
    tmap = assignTextons(fbRun(textonData.fb,im),textonData.tex);
    [tg,gtheta] = tgmo(tmap,k,idiag*radius(4),norient,...
                      'smooth','savgol','sigma',idiag*radius(4));
                  
    % nonmax suppression and max over orientations
    [tg_max,maxo] = max(tg,[],3);
    pb = zeros(h,w);
    theta = zeros(h,w);
    r = 2.5;
    for j = 1:norient,
      mask = (maxo == j);
      a = fitparab(tg(:,:,j),r,r,gtheta(j));
      pbi = nonmax(max(0,a),gtheta(j));
      pb = max(pb,pbi.*mask);
      theta = theta.*~mask + gtheta(j).*mask;
    end
    pb = max(0,min(1,pb));

    % mask out 1-pixel border where nonmax suppression fails
    pb(1,:) = 0;
    pb(end,:) = 0;
    pb(:,1) = 0;
    pb(:,end) = 0;              
%     imshow(pb, 'border', 'tight');
%     imwrite(tg_max, [dst_path img_files(i).name(1:end-4) '_tg_response.png'], 'png')
%     imwrite(pb, [dst_path img_files(i).name(1:end-4) '_tg_bry.png'], 'png')
    
    

% Use NMS to detect the edgels - this NMS routine ouputs subpixel edgel tokens
margin = 4;
dirx = cos(mod(theta+pi/2, pi));
diry = sin(mod(theta+pi/2, pi));
[subpix_x, subpix_y, subpix_dir_x, subpix_dir_y, subpix_grad_x, subpix_grad_y] = NMS_token(dirx, diry,pb, pb>thresh, margin);

% magnitude of the gradient at the maxima
mag_e = sqrt(subpix_grad_x.^2+subpix_grad_y.^2);

Hx = pb.*dirx;
Hy = pb.*diry;
sigma = 0.7;
Px   = resampled_filter_2d(pb, @Gx_2d_op, sigma, 0);
Py   = resampled_filter_2d(pb, @Gy_2d_op, sigma, 0);
Pxx  = resampled_filter_2d(pb, @Gxx_2d_op, sigma, 0);
Pxy  = resampled_filter_2d(pb, @Gxy_2d_op, sigma, 0);
Pyy  = resampled_filter_2d(pb, @Gyy_2d_op, sigma, 0);

Hxx  = resampled_filter_2d(Hx, @Gx_2d_op, sigma, 0);
Hyy  = resampled_filter_2d(Hy, @Gy_2d_op, sigma, 0);
Hxy  = resampled_filter_2d(Hy, @Gx_2d_op, sigma, 0);

[yd, xd] = size(Hx);
[xx,yy] = meshgrid(1:xd,1:yd);
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

%normalize
F_mag = sqrt(Fx_e.^2 + Fy_e.^2);
Fx_e = Fx_e./F_mag;
Fy_e = Fy_e./F_mag;

% Edge Direction (tangent to the level set) is orthogonal to the gradient
subpix_dir_x2 = (-Fy_e);
subpix_dir_y2 = (Fx_e);

% construct the edge map
edg  = [subpix_x' subpix_y' atan2(subpix_dir_y2', subpix_dir_x2') mag_e'];