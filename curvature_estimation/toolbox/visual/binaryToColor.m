function [ colorMap ] = binaryToColor( binaryMap, color, background )
% Visualize a binary image on a background image or visualize multiple 
% binary images with different color by calling this function several
% times.
%
% USAGE
%  [ colorMap ] = binaryToColor( binaryMap, color )
%  [ colorMap ] = binaryToColor( binaryMap, color, background )
%
% INPUTS
%  binaryMap  - [h x w] input binary image
%  color      - [char/1 x 3] the color used to draw the binary image
%  background - [h x w x 3] background image
%
% OUTPUTS
%  colorMap   - [h x w] viualization results
%
% Code written by Xiaoyan Li, 2014.

if nargin<3
    background = ones([size(binaryMap) 3]);
end
% reture gray image in specified color
if ischar(color)
    key = {'y','m','c','r','g','b','k','o','w'};
    value = {[1 1 0],[1 0 1],[0 1 1],[1 0 0],[0 1 0],[0 0 1],[0 0 0], [1 .5 0], [1 1 1]};
    colorMap = containers.Map(key,value);
    if ~isKey(colorMap,color), error('this color name is not supported.'); end
    colorArray = colorMap(color);
else
    if (numel(color) == 3) && all(color<=1)
        colorArray = color;
    else
        error('Proper color values are not provided.')
    end
end

colorMap = im2double(background);
bin = logical(binaryMap);
for i = 1:3
%     if ~islogical(binaryMap)
        addimg = (1-binaryMap*(1-colorArray(i))).*bin;
%     else
%         addimg = bin*colorArray(i);
%     end
    colorMap(:,:,i) = (~bin).*double(colorMap(:,:,i))+addimg;
end
end
