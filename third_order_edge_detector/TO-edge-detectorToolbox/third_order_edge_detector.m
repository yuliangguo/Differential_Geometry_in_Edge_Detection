function [TO_edge_map, gen_edge_map, grad_mag] = third_order_edge_detector(img, sigma, n, threshold, QaD_flag)
%function [TO_edge_map, gen_edge_map] = third_order_edge_detector(img, sigma, n, threshold,QaD_flag)
% 
% --------------------------------
% Third-order edge detector (v1.0)
% --------------------------------
%
% This work first appeared in :
%   Tamrakar, A. and Kimia, B. B., "No grouping left behind: From edges to
%   curve fragments", ICCV, 2007
%
% Summary:
%   First compute the gradient image and perform NMS to locate edges.
%   Next locate the dubpixel location of the maximum and compute the
%   orientation of the edge curve using third-order derivatives.
%
%   F(x,y) = grad(I)^T * Hessian(I) * grad(I) = 0; (defines the edge curve)
%
%   So the edge orientation needs to be 
%       [Fx, Fy] = grad(F)/|grad(F)|
%   evaluated at the zero level set of F(x,y).
%
% Parameters:
%   img:       input image
%   sigma:     scale of the operator (>0.7), typically 1
%   n:         interpolating factor (2^n number of subdivions within the
%   pixel) typically 1
%   threshold: gradient magnitude threshold (not normalized, typically 2)
%   QaD_flag:  Quck-and-Dirty flag: compute subpixel derivatives by
%              interpolation instead of resampled fitlers
% Output
%   TO_edge_map:  Third-order edge orientation map
%   gen_edge_map: Generic Edge map
%
%       Note: An edgemap is a matrix of detected edgels [each row contains [x, y, theta,
%       strength]
%
% See also resampled_filt_2d
%
% (c) LEMS, Brown University
% Amir Tamrakar (amir_tamrakar@brown.edu)
% October 2007

%% Preprocess
if (size(img,3)>1) %convert from rgb to gray
    img = rgb2gray(img);
end
img = double(img);
margin = ceil(4*sigma);

%% pad the image
[h,w] = size(img); 
img = padarray(img,[margin margin],'symmetric','both');

%% Compute derivatives using Shifted Gaussian filters
Ix   = resampled_filter_2d(img, @Gx_2d_op, sigma, n);
Iy   = resampled_filter_2d(img, @Gy_2d_op, sigma, n);

%% Compute supixel position
grad_mag = sqrt(Ix.^2+Iy.^2); %gradient magnitude

% compute F(x,y)
% F = Ix.^2.*Ixx + 2*Ix.*Iy.*Ixy + Iy.^2.*Iyy;

% Note: we can detect the edges either by finding maxima of grad_mag or
% zero crossings of F

% We'll use NMS to detect the edgels - this NMS routine ouputs subpixel edgel tokens
[subpix_x, subpix_y, subpix_dir_x, subpix_dir_y, subpix_grad_x, subpix_grad_y] = NMS_token(Ix, Iy, grad_mag, grad_mag>threshold, margin);

%% Compute the corrected orientation
%
%  The orientations of the edges computed from the gradient are incorrect. 
%  Nevertheless the subpixel locations computed using these orientations
%  should be decent. So using the subpixel locations, we can recompute the 
%  orientations using the third-order formulation.

if (QaD_flag)
    % Quick and dirty way of computing the derivatives at the edges
    %  --> interpolate the derivatives at the edge locations
    
    Ixx  = resampled_filter_2d(img, @Gxx_2d_op, sigma, n);
    Iyy  = resampled_filter_2d(img, @Gyy_2d_op, sigma, n);
    Ixy  = resampled_filter_2d(img, @Gxy_2d_op, sigma, n);
    Ixxy = resampled_filter_2d(img, @Gxxy_2d_op, sigma, n);
    Ixyy = resampled_filter_2d(img, @Gxyy_2d_op, sigma, n);
    Ixxx = resampled_filter_2d(img, @Gxxx_2d_op, sigma, n);
    Iyyy = resampled_filter_2d(img, @Gyyy_2d_op, sigma, n);

    [yd, xd] = size(Ix);
    [xx,yy] = meshgrid(1:xd,1:yd);
    Ix_e   = interp2(xx, yy, Ix,   subpix_x, subpix_y, 'cubic');
    Iy_e   = interp2(xx, yy, Iy,   subpix_x, subpix_y, 'cubic');
    Ixx_e  = interp2(xx, yy, Ixx,  subpix_x, subpix_y, 'cubic');
    Iyy_e  = interp2(xx, yy, Iyy,  subpix_x, subpix_y, 'cubic');
    Ixy_e  = interp2(xx, yy, Ixy,  subpix_x, subpix_y, 'cubic');
    Ixxy_e = interp2(xx, yy, Ixxy, subpix_x, subpix_y, 'cubic');
    Ixyy_e = interp2(xx, yy, Ixyy, subpix_x, subpix_y, 'cubic');
    Ixxx_e = interp2(xx, yy, Ixxx, subpix_x, subpix_y, 'cubic');
    Iyyy_e = interp2(xx, yy, Iyyy, subpix_x, subpix_y, 'cubic');
else
    % compute the derivatives by resampling the fitlers at the zero
    % crossings
    
    % I'll update this soon. Keep using the QaD flag for now.
    
end

% magnitude of the gradient at the maxima
mag_e = sqrt(subpix_grad_x.^2+subpix_grad_y.^2);

% compute [Fx, Fy] at zero crossings
Fx_e = 2*Ix_e.*Ixx_e.^2 + 2*Ix_e.*Ixy_e.^2 + 2*Iy_e.*Ixx_e.*Ixy_e + 2*Iy_e.*Iyy_e.*Ixy_e + 2*Ix_e.*Iy_e.*Ixxy_e + Iy_e.^2.*Ixyy_e + Ix_e.^2.*Ixxx_e;
Fy_e = 2*Iy_e.*Iyy_e.^2 + 2*Iy_e.*Ixy_e.^2 + 2*Ix_e.*Ixx_e.*Ixy_e + 2*Ix_e.*Iyy_e.*Ixy_e + 2*Ix_e.*Iy_e.*Ixyy_e + Ix_e.^2.*Ixxy_e + Iy_e.^2.*Iyyy_e;

%normalize
F_mag = sqrt(Fx_e.^2 + Fy_e.^2);
Fx_e = Fx_e./F_mag;
Fy_e = Fy_e./F_mag;

% Edge Direction (tangent to the level set) is orthogonal to the gradient
subpix_dir_x2 = (-Fy_e);
subpix_dir_y2 = (Fx_e);

amp = 2^n;
% construct the edge map
gen_edge_map = [(subpix_x'-1)/amp (subpix_y'-1)/amp atan2(subpix_dir_y',  subpix_dir_x')  mag_e'];
TO_edge_map  = [(subpix_x'-1)/amp (subpix_y'-1)/amp atan2(subpix_dir_y2', subpix_dir_x2') mag_e'];
ind = (TO_edge_map(:,1)>margin-0.5) & (TO_edge_map(:,1)<margin-0.5+w) & ...
    (TO_edge_map(:,2)>margin-0.5) & (TO_edge_map(:,2)<margin-0.5+h);
gen_edge_map = gen_edge_map(ind,:); gen_edge_map(:,1:2) = gen_edge_map(:,1:2)-margin+0.5;
TO_edge_map = TO_edge_map(ind,:); TO_edge_map(:,1:2) = TO_edge_map(:,1:2)-margin+0.5;

%% Display edgels (for debug only )
% figure;
% 
% axis([120 180 120 180]);
% % imshow(img, [0 255]);
% % display edges computed by the gradient operator
% %disp_edg(gen_edge_map, 'r');
% % display edges whose orientations have been computed by the third-order operator
% disp_edg(TO_edge_map, 'g');
grad_mag = imresize(grad_mag,size(img)+2*margin);
grad_mag = grad_mag(margin+1:h+margin,margin+1:w+margin);