function binary_map = create_shape_on_binary_map( fxy, background )
[h,w] = size(background);
[xx,yy] = meshgrid(1:w,1:h);
Fxy = fxy(xx,yy);
background(Fxy<0) = 1;
binary_map = background;
end

