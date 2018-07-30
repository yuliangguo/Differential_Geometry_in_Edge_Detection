function visualize_edge_matches( handle, edgemap1, edgemap2, thetamap1, thetamap2, match, box, img, thresh )
[h,w] = size(match);
ind1 = find(match);
ind2 = match(ind1);
match2_bin = zeros(h,w);
match2_bin(ind2) = 1;

if nargin<9
    thresh = 0;
end
theta1 = wrapToPi(thetamap1(ind1));
theta2 = wrapToPi(thetamap2(ind2));
dif = theta1-theta2;
thetadiff = min(cat(3,abs(dif),abs(dif+pi),abs(dif-pi)),[],3);
mask1 = zeros(h,w);
mask1(ind1(thetadiff>=thresh)) = 1;
mask2 = zeros(h,w);
mask2(ind2(thetadiff>=thresh)) = 1;

mask = zeros(h,w);
mask(box(3):box(4),box(1):box(2)) = 1;
edgemap1 = edgemap1.*mask.*mask1;
edgemap2 = edgemap2.*mask.*mask2;

[ edg1 ] = convert_to_edg( edgemap1, thetamap1 );
if nargin>7
   disp_edg_in_box(handle, edg1, box, 'g', img);
else
    disp_edg_in_box(handle, edg1, box, 'g');
end
[ edg2 ] = convert_to_edg( edgemap2, thetamap2 );
disp_edg_in_box(handle, edg2, box, 'r');
    
ind1 = find(match.*mask.*mask1);
ind2 = match(ind1);
[i1,j1] = ind2sub([h,w],ind1);
[i2,j2] = ind2sub([h,w],ind2);
figure(handle);

axis(box);
hold on;
for e = 1:length(i1)
    plot([j1(e) j2(e)],[i1(e) i2(e)],'b');
end
hold off;
% axis equal;
axis off;


end

