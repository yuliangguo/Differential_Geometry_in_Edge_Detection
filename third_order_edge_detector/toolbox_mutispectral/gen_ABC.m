function [A, Ax, Ay, Axx, Axy, Ayy, B, Bx, By, Bxx, Bxy, Byy, C, Cx, Cy, Cxx, Cxy, Cyy] = gen_ABC( img, n, sigma )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
[ Ix1, Iy1, Ixx1, Ixy1, Iyy1, Ixxx1, Ixxy1, Ixyy1, Iyyy1 ] = gen_derivative( img(:,:,1), n, sigma );
[ Ix2, Iy2, Ixx2, Ixy2, Iyy2, Ixxx2, Ixxy2, Ixyy2, Iyyy2 ] = gen_derivative( img(:,:,2), n, sigma );
[ Ix3, Iy3, Ixx3, Ixy3, Iyy3, Ixxx3, Ixxy3, Ixyy3, Iyyy3 ] = gen_derivative( img(:,:,3), n, sigma );

A   = Ix1.^2 + Ix2.^2 + Ix3.^2;

Ax  = 2*Ix1.*Ixx1 + 2*Ix2.*Ixx2 + 2*Ix3.*Ixx3;

Ay  = 2*Ix1.*Ixy1 + 2*Ix2.*Ixy2 + 2*Ix3.*Ixy3;

Axx = 2*Ix1.*Ixxx1 + 2*Ixx1.^2 + 2*Ix2.*Ixxx2 + 2*Ixx2.^2 + 2*Ix3.*Ixxx3 + 2*Ixx3.^2;

Axy = 2*Ix1.*Ixxy1 + 2*Ixy1.*Ixx1 + 2*Ix2.*Ixxy2 + 2*Ixy2.*Ixx2 + 2*Ix3.*Ixxy3 + 2*Ixy3.*Ixx3;

Ayy = 2*Ix1.*Ixyy1 + 2*Ixy1.^2 + 2*Ix2.*Ixyy2 + 2*Ixy2.^2 + 2*Ix3.*Ixyy3 + 2*Ixy3.^2;

B   = Ix1.*Iy1+Ix2.*Iy2+Ix3.*Iy3;

Bx  = Ix1.*Ixy1 + Iy1.*Ixx1 + Ix2.*Ixy2 + Iy2.*Ixx2 + Ix3.*Ixy3 + Iy3.*Ixx3;

By  = Ix1.*Iyy1 + Iy1.*Ixy1 + Ix2.*Iyy2 + Iy2.*Ixy2 + Ix3.*Iyy3 + Iy3.*Ixy3;

Bxx = Ix1.*Ixxy1 + Ixx1.*Ixy1 + Iy1.*Ixxx1 + Ixy1.*Ixx1 + Ix2.*Ixxy2 + Ixx2.*Ixy2 ...
     + Iy2.*Ixxx2 + Ixy2.*Ixx2 + Ix3.*Ixxy3 + Ixx3.*Ixy3 + Iy3.*Ixxx3 + Ixy3.*Ixx3;
 
Bxy = Ix1.*Ixyy1 + Ixy1.^2 + Iy1.*Ixxy1 + Iyy1.*Ixx1 + Ix2.*Ixyy2 + Ixy2.^2 ...
     + Iy2.*Ixxy2 + Iyy2.*Ixx2 + Ix3.*Ixyy3 + Ixy3.^2 + Iy3.*Ixxy3 + Iyy3.*Ixx3;
 
Byy = Ix1.*Iyyy1 + Ixy1.*Iyy1 + Iy1.*Ixyy1 + Ix2.*Iyyy2 + Ixy2.*Iyy2 ...
     + Iy2.*Ixyy2 + Ix3.*Iyyy3 + Ixy3.*Iyy3 + Iy3.*Ixyy3;
 
C   = Iy1.^2+Iy2.^2+Iy3.^2;

Cx  = 2*Iy1.*Ixy1 + 2*Iy2.*Ixy2 + 2*Iy3.*Ixy3;

Cy  = 2*Iy1.*Iyy1 + 2*Iy2.*Iyy2 + 2*Iy3.*Iyy3;

Cxx = 2*Iy1.*Ixxy1 + 2*Ixy1.^2 + 2*Iy2.*Ixxy2 + 2*Ixy2.^2 + 2*Iy3.*Ixxy3 + 2*Ixy3.^2;

Cxy = 2*Iy1.*Ixyy1 + 2*Iyy1.*Ixy1 + 2*Iy2.*Ixyy2 + 2*Iyy2.*Ixy2 + 2*Iy3.*Ixyy3 + 2*Iyy3.*Ixy3;

Cyy = 2*Iy1.*Iyyy1 + 2*Iyy1.^2 + 2*Iy2.*Iyyy2 + 2*Iyy2.^2 + 2*Iy3.*Iyyy3 + 2*Iyy3.^2;

end

