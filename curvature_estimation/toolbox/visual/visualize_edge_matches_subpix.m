function visualize_edge_matches_subpix( handle, edgeinfo1, edgeinfo2, match, box, img )
assert(size(match,2) == 2);
pos1 = edgeinfo1(match(:,1),1:4);
pos2 = edgeinfo2(match(:,2),1:4);
ind = find(pos1(:,1)>box(1) & pos1(:,1)<box(2) & pos1(:,2)>box(3) & pos1(:,2)<box(4) & ...
    pos2(:,1)>box(1) & pos2(:,1)<box(2) & pos2(:,2)>box(3) & pos2(:,2)<box(4));
pos1 = pos1(ind,:);
pos2 = pos2(ind,:);

figure(handle);
if nargin>5
    imshow(img,'border','tight');
end
axis(box);
disp_edg(pos1,'g');
disp_edg(pos2,'r');
hold on;
for e = 1:size(pos1,1)
    plot([pos1(e,1)+1 pos2(e,1)+1],[pos1(e,2)+1 pos2(e,2)+1],'b');
end
hold off;
end

