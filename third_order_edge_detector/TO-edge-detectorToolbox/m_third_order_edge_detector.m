function [TO_edge_map, gen_edge_map] = m_third_order_edge_detector(img, weights, sigma, threshold)
% input:
% img:          original image, only deal with gray map now
% weights:      weights for multi-scale maps
% sigma:        a vector indicating multi-scale gaussian filters. This sigama should be the same as used for training weights
% n:            2^n number of subdivions within the pixel
% threshold:    gradient threshold for selecting output edges
% QaD_flag:     Quick and dirty way of computing the derivatives at the edges

% now just fix these two param
n = 1;
QaD_flag = 1;

if (size(img,3)>1) %convert from rgb to gray
    img = rgb2gray(img);
end
img = double(img);
    
% compute gaussian derivation response map for each scale
    for s = 1:length(sigma)
        Ix(:,:,s)   = resampled_filter_2d(img, @Gx_2d_op, sigma(s), n);
        Iy(:,:,s)   = resampled_filter_2d(img, @Gy_2d_op, sigma(s), n);
        grad_mag(:,:,s) = sqrt(Ix(:,:,s).^2+Iy(:,:,s).^2);
    end
    [h, w, ~] = size(Ix);

    new_weights = repmat(weights, [1, h, w]);
    new_weights = shiftdim(new_weights, 1);

    % combine multiscale dirivation map
    Ix = Ix.*new_weights;
    Iy = Iy.*new_weights;
    Ix = sum(Ix, 3);
    Iy = sum(Iy, 3);

    %% Compute 
    % % perform the same normalization
    % for h1 = 1:size(grad_mag,1)
    %     for w1 = 1:size(grad_mag,2)
    %         grad_mag(h1,w1,:) = grad_mag(h1,w1,:)/norm(shiftdim(grad_mag(h1,w1,:),2));
    %     end
    % end

    grad_mag = grad_mag.*new_weights; %gradient magnitude
    grad_mag = sum(grad_mag,3);

    % compute F(x,y)
    % F = Ix.^2.*Ixx + 2*Ix.*Iy.*Ixy + Iy.^2.*Iyy;

    % Note: we can detect the edges either by finding maxima of grad_mag or
    % zero crossings of F

    % We'll use NMS to detect the edgels - this NMS routine ouputs subpixel edgel tokens
    margin = ceil(4*sigma(1));
    [subpix_x, subpix_y, subpix_dir_x, subpix_dir_y, subpix_grad_x, subpix_grad_y] = NMS_token(Ix, Iy, grad_mag, grad_mag>threshold, margin);

    % magnitude of the gradient at the maxima
    mag_e = sqrt(subpix_grad_x.^2+subpix_grad_y.^2);
    gen_edge_map = [subpix_x'/(n+1) subpix_y'/(n+1) atan2(subpix_dir_y',  subpix_dir_x')  mag_e'];

    %% Compute the corrected orientation
    %
    %  The orientations of the edges computed from the gradient are incorrect. 
    %  Nevertheless the subpixel locations computed using these orientations
    %  should be decent. So using the subpixel locations, we can recompute the 
    %  orientations using the third-order formulation.

    if (QaD_flag)
        % Quick and dirty way of computing the derivatives at the edges
        %  --> interpolate the derivatives at the edge locations
        for s = 1:length(sigma)
            Ixx(:,:,s)  = resampled_filter_2d(img, @Gxx_2d_op, sigma(s), n);
            Iyy(:,:,s)  = resampled_filter_2d(img, @Gyy_2d_op, sigma(s), n);
            Ixy(:,:,s)  = resampled_filter_2d(img, @Gxy_2d_op, sigma(s), n);
            Ixxy(:,:,s) = resampled_filter_2d(img, @Gxxy_2d_op, sigma(s), n);
            Ixyy(:,:,s) = resampled_filter_2d(img, @Gxyy_2d_op, sigma(s), n);
            Ixxx(:,:,s) = resampled_filter_2d(img, @Gxxx_2d_op, sigma(s), n);
            Iyyy(:,:,s) = resampled_filter_2d(img, @Gyyy_2d_op, sigma(s), n);
        end

        % combine with the same weights
        Ixx = sum(Ixx.*new_weights, 3);
        Iyy = sum(Iyy.*new_weights, 3);
        Ixy = sum(Ixy.*new_weights, 3);
        Ixxy = sum(Ixxy.*new_weights, 3);
        Ixyy = sum(Ixyy.*new_weights, 3);
        Ixxx = sum(Ixxx.*new_weights, 3);
        Iyyy = sum(Iyyy.*new_weights, 3);

        % interp
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

    % construct the edge map
    TO_edge_map  = [subpix_x'/(n+1) subpix_y'/(n+1) atan2(subpix_dir_y2', subpix_dir_x2') mag_e'];

    % %% Display edgels (for debug only )
    % figure;
    % imshow(img, [0 255]);
    % % display edges computed by the gradient operator
    % disp_edg(gen_edge_map, 'r');
    % % % display edges whose orientations have been computed by the third-order operator
    hold on;
    axis([120 180 120 180]);
    disp_edg(TO_edge_map, 'g');
    hold off;
    end
