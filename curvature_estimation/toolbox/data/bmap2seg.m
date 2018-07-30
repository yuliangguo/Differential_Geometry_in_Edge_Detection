function [seg] = bmap2seg(bmap,width,height)

if nargin<3,
  [height,width] = size(bmap);
end
[h,w] = size(bmap);

filled_bmap = imfill(bmap,'holes');
filled_bmap(bmap) = 0;

[s] = bwlabel(filled_bmap);

% check width and height
ar1 = width / height;
ar2 = w / h;
if width>w || height>h || abs(ar1-ar2)>0.01,
  error(sprintf('Can''t convert %dx%d bmap to %dx%d seg.',w,h,width,height));
end

  
if w==width && h==height,
  
  seg = s;

else
  
  seg = zeros(height,width);
  for x = 1:w,
    for y = 1:h,
      if s(y,x)>0,
        j = 1+floor((y-1)*height/h);
        i = 1+floor((x-1)*width/w);
        bmap(j,i) = s(y,x);
      end
    end
  end

end

