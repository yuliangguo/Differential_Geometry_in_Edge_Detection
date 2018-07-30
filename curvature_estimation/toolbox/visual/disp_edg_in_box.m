function disp_edg_in_box( h, edginfo, box, col, img )
figure(h);
if nargin>4
    imshow(img, 'border', 'tight');
end
axis(box);
% note box use matlab coordinates, while edginfo use c++ coordinates
ind = find(edginfo(:,1)>box(1)-1 & edginfo(:,1)<box(2)-1 & edginfo(:,2)>box(3)-1 ...
    & edginfo(:,2)<box(4)-1);
disp_edg(edginfo(ind,:),col)

end

