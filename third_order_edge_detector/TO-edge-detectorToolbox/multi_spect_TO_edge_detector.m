function [TO_edge, gen_edge]=multi_spect_TO_edge_detector(img, interp, sigma, T0, col_space)
%% preprocess: change the image into specified color space 
if strcmp(col_space, 'Lab')
    img_n = rgb2lab(img);
elseif strcmp(col_space, 'HSV')
    img_n = rgb2hsv(img);
elseif strcmp(col_space, 'oppo')
    img_p = zeros(size(img));
    img0 = double(img);
    img_p(:,:,3)=(img0(:,:,1)+img0(:,:,2)+img0(:,:,3))/sqrt(3);
    img_p(:,:,1)=(img0(:,:,1)-img0(:,:,2))/sqrt(2);
    img_p(:,:,2)=(img0(:,:,1)+img0(:,:,2)-2*img0(:,:,3))/sqrt(6);
    img_n = img_p;
else
    img_n= double(img);
end
margin = ceil(4*sigma);

%% pad the image
[h,w,~] = size(img); 
img_n = padarray(img_n,[margin margin,0],'symmetric','both');

%% compute derivatives using Shifted Gaussian filters
[A, Ax, Ay, Axx, Axy, Ayy, B, Bx, By, Bxx, Bxy, Byy, C, Cx, Cy, Cxx, Cxy, Cyy] = gen_ABC( img_n, interp, sigma );

%% Compute supixel position
Lambada = (A+C+sqrt((A-C).^2+4*B.^2))/2;
ind1 = find(B>10^-2);
ind2 = find(B<=10^-2);
Eta1(ind1) = B(ind1)./sqrt(B(ind1).^2+(Lambada(ind1)-A(ind1)).^2);
Eta2(ind1) = (Lambada(ind1)-A(ind1))./sqrt(B(ind1).^2+(Lambada(ind1)-A(ind1)).^2);
Eta1(ind2) = (Lambada(ind2)-C(ind2))./sqrt(B(ind2).^2+(Lambada(ind2)-C(ind2)).^2);
Eta2(ind2) = B(ind2)./sqrt(B(ind2).^2+(Lambada(ind2)-C(ind2)).^2);
Eta1 = reshape(Eta1,size(Lambada));
Eta2 = reshape(Eta2,size(Lambada));

% maxLam = max(max(Lambada));
% minLam = min(min(Lambada));
Lambada = sqrt(Lambada/3);
[subpix_x, subpix_y, subpix_dir_x, subpix_dir_y, subpix_grad_x, subpix_grad_y] = NMS_token(Eta1, Eta2, ...
    Lambada, Lambada>T0, margin);

%% Compute the corrected orientation
[yd, xd] = size(A);
[xx,yy] = meshgrid(1:xd,1:yd);
A_m    = interp2(xx, yy, A,    subpix_x, subpix_y, 'cubic');
Ax_m   = interp2(xx, yy, Ax,   subpix_x, subpix_y, 'cubic');
Ay_m   = interp2(xx, yy, Ay,   subpix_x, subpix_y, 'cubic');
Axx_m  = interp2(xx, yy, Axx,  subpix_x, subpix_y, 'cubic');
Ayy_m  = interp2(xx, yy, Ayy,  subpix_x, subpix_y, 'cubic');
Axy_m  = interp2(xx, yy, Axy,  subpix_x, subpix_y, 'cubic');
B_m    = interp2(xx, yy, B,    subpix_x, subpix_y, 'cubic');
Bx_m   = interp2(xx, yy, Bx,   subpix_x, subpix_y, 'cubic');
By_m   = interp2(xx, yy, By,   subpix_x, subpix_y, 'cubic');
Bxx_m  = interp2(xx, yy, Bxx,  subpix_x, subpix_y, 'cubic');
Byy_m  = interp2(xx, yy, Byy,  subpix_x, subpix_y, 'cubic');
Bxy_m  = interp2(xx, yy, Bxy,  subpix_x, subpix_y, 'cubic');
C_m    = interp2(xx, yy, C,    subpix_x, subpix_y, 'cubic');
Cx_m   = interp2(xx, yy, Cx,   subpix_x, subpix_y, 'cubic');
Cy_m   = interp2(xx, yy, Cy,   subpix_x, subpix_y, 'cubic');
Cxx_m  = interp2(xx, yy, Cxx,  subpix_x, subpix_y, 'cubic');
Cyy_m  = interp2(xx, yy, Cyy,  subpix_x, subpix_y, 'cubic');
Cxy_m  = interp2(xx, yy, Cxy,  subpix_x, subpix_y, 'cubic');
Eta1_m = interp2(xx, yy, Eta1, subpix_x, subpix_y, 'cubic');
Eta2_m = interp2(xx, yy, Eta2, subpix_x, subpix_y, 'cubic');

Eta1x_m = ((Ax_m-Cx_m).*((A_m-C_m).^2+4*B_m.^2) - (A_m-C_m).*((A_m-C_m).*(Ax_m-Cx_m)+4.*B_m.*Bx_m)) ...
          ./(4*Eta1_m.*((A_m-C_m).^2+4*B_m.^2).^1.5);
Eta2x_m = -((Ax_m-Cx_m).*((A_m-C_m).^2+4*B_m.^2) - (A_m-C_m).*((A_m-C_m).*(Ax_m-Cx_m)+4.*B_m.*Bx_m)) ...
          ./(4*Eta2_m.*((A_m-C_m).^2+4*B_m.^2).^1.5);
Eta1y_m = ((Ay_m-Cy_m).*((A_m-C_m).^2+4*B_m.^2) - (A_m-C_m).*((A_m-C_m).*(Ay_m-Cy_m)+4.*B_m.*By_m)) ...
          ./(4*Eta1_m.*((A_m-C_m).^2+4*B_m.^2).^1.5);
Eta2y_m = -((Ay_m-Cy_m).*((A_m-C_m).^2+4*B_m.^2) - (A_m-C_m).*((A_m-C_m).*(Ay_m-Cy_m)+4.*B_m.*By_m)) ...
          ./(4*Eta2_m.*((A_m-C_m).^2+4*B_m.^2).^1.5);

Fx = Axx_m.*Eta1_m.^3 + 3*Ax_m.*Eta1_m.^2.*Eta1x_m + Eta1_m.^2.*Eta2_m.*Axy_m + 2*Eta1_m.^2.*Eta2_m.*Bxx_m ...
    + 2*Eta1_m.*Eta2_m.*Eta1x_m.*Ay_m + 4*Eta1_m.*Eta2_m.*Eta1x_m.*Bx_m + Eta1_m.^2.*Eta2x_m.*Ay_m + ...
    2*Eta1_m.^2.*Eta2x_m.*Bx_m + Eta1_m.*Eta2_m.^2.*Cxx_m + 2*Eta1_m.*Eta2_m.^2.*Bxy_m + ...
    Eta1x_m.*Eta2_m.^2.*Cx_m + 2*Eta1x_m.*Eta2_m.^2.*By_m + 2*Eta1_m.*Eta2_m.*Eta2x_m.*Cx_m + ...
    4*Eta1_m.*Eta2_m.*Eta2x_m.*By_m + Cxy_m.*Eta2_m.^3 + 3*Cy_m.*Eta2_m.^2.*Eta2x_m;

Fy = Axy_m.*Eta1_m.^3 + 3*Ax_m.*Eta1_m.^2.*Eta1y_m + Eta1_m.^2.*Eta2_m.*Ayy_m + 2*Eta1_m.^2.*Eta2_m.*Bxy_m ...
    + 2*Eta1_m.*Eta2_m.*Eta1y_m.*Ay_m + 4*Eta1_m.*Eta2_m.*Eta1y_m.*Bx_m + Eta1_m.^2.*Eta2y_m.*Ay_m + ...
    2*Eta1_m.^2.*Eta2y_m.*Bx_m + Eta1_m.*Eta2_m.^2.*Cxy_m + 2*Eta1_m.*Eta2_m.^2.*Byy_m + ...
    Eta1y_m.*Eta2_m.^2.*Cx_m + 2*Eta1y_m.*Eta2_m.^2.*By_m + 2*Eta1_m.*Eta2_m.*Eta2y_m.*Cx_m + ...
    4*Eta1_m.*Eta2_m.*Eta2y_m.*By_m + Cyy_m.*Eta2_m.^3 + 3*Cy_m.*Eta2_m.^2.*Eta2y_m;

F_mag = sqrt(Fx.^2 + Fy.^2);
Fx = Fx./F_mag;
Fy = Fy./F_mag;

dir_x = (-Fy);
dir_y = (Fx);

angle = atan2(dir_y, dir_x);
% angle0 = atan2(subpix_dir_y, subpix_dir_x);

% figure;
% imshow(img);
% axis([120 180 120 180]);
amp=2^interp;
Mag_e = sqrt(subpix_grad_x.^2+subpix_grad_y.^2);
% ind = sub2ind(size(Lambada),round(subpix_y),round(subpix_x));
TO_edge  = [(subpix_x'-1)/amp (subpix_y'-1)/amp angle' Mag_e']; % c++ coordinate
gen_edge = [(subpix_x'-1)/amp (subpix_y'-1)/amp atan2(subpix_dir_y,subpix_dir_x)' Mag_e'];
ind = (TO_edge(:,1)>margin-0.5) & (TO_edge(:,1)<margin-0.5+w) & ...
    (TO_edge(:,2)>margin-0.5) & (TO_edge(:,2)<margin-0.5+h);
gen_edge = gen_edge(ind,:); gen_edge(:,1:2) = gen_edge(:,1:2)-margin;
TO_edge = TO_edge(ind,:); TO_edge(:,1:2) = TO_edge(:,1:2)-margin;

end
