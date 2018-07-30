function [ Gyyy ] = Gaussian_yyy( sigma, dx,dy )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
half_size = ceil(sigma*4+max(abs(dx),abs(dy)));
[X,Y] = meshgrid(-half_size:half_size,-half_size:half_size);
Gyyy = (3*(Y-dy)*sigma^2-(Y-dy).^3).*exp(-((X-dx).^2+(Y-dy).^2)/(2*sigma^2))/(2*pi*sigma^8); 
end