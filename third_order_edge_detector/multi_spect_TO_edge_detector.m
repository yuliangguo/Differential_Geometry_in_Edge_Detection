function [TO_edge]=multi_spect_TO_edge_detector(img, interp, sigma, TH, col_space, is_vis)
%figure,imshow(img);
img_n = zeros(size(img));
if strcmp(col_space, 'Lab')
    colorTransform = makecform('srgb2lab');
    img_p = applycform(img,colorTransform);
    img_p = double(img_p);
%     img_p = RGB2Lab(img);
    img_n(:,:,1)=img_p(:,:,1)/100;
    img_n(:,:,2)=(img_p(:,:,2)+128)/255;
    img_n(:,:,3)=(img_p(:,:,3)+128)/255;
elseif strcmp(col_space, 'HSV')
    img_p = rgb2hsv(img);
    img_n = double(img_p);
elseif strcmp(col_space, 'oppo')
    img_p = zeros(size(img));
    img0 = double(img);
    img_p(:,:,3)=(img0(:,:,1)+img0(:,:,2)+img0(:,:,3))/sqrt(3);
    img_p(:,:,1)=(img0(:,:,1)-img0(:,:,2))/sqrt(2);
    img_p(:,:,2)=(img0(:,:,1)+img0(:,:,2)-2*img0(:,:,3))/sqrt(6);
    img_n(:,:,3)=img_p(:,:,3)/sqrt(3)/255;
    img_n(:,:,1)=(img_p(:,:,1)*sqrt(2)+255)/2/255;
    img_n(:,:,2)=(img_p(:,:,2)*sqrt(6)+2*255)/4/255;
else
    img_n= im2double(img);
end
% [y,x,z]=size(img);
% n=1; %interpolation 2^n
% sigma=0.7;
% T0=0.01;
n=interp;

[A, Ax, Ay, Axx, Axy, Ayy, B, Bx, By, Bxx, Bxy, Byy, C, Cx, Cy, Cxx, Cxy, Cyy] = gen_ABC( img_n, n, sigma );

Lambada = (A+C+sqrt((A-C).^2+4*B.^2))/2;
Eta1    = sqrt((1+(A-C)./sqrt((A-C).^2+4*B.^2))/2);
Eta2    = sign(B).*sqrt((1-(A-C)./sqrt((A-C).^2+4*B.^2))/2);

maxLam = max(max(Lambada));
minLam = min(min(Lambada));
margin = ceil(4*sigma);
Lambada0 = (Lambada-minLam)*40/(maxLam-minLam);
[subpix_x, subpix_y, subpix_dir_x, subpix_dir_y, subpix_grad_x, subpix_grad_y] = NMS_token(Eta1, Eta2, ...
    Lambada0, Lambada0>TH, margin);

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

if(is_vis)
figure;
imshow(img, 'border', 'tight');
end
% axis([120 180 120 180]);
length=1;
amp=2^n;
% edge_map( subpix_x/amp, subpix_y/amp, angle0, length, 'r' );
if(is_vis)
edge_map( subpix_x/amp, subpix_y/amp, angle,  length, 'g' );
end
% edgex = subpix_x/amp;
% edgey = subpix_x/amp;
% orientation = angle;
Mag_e = sqrt(subpix_grad_x.^2+subpix_grad_y.^2);
TO_edge  = [subpix_x'/amp subpix_y'/amp angle' Mag_e'];
end
