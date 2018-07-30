function [orientation_map, edge_map]= gt_orientation(human_seg)

orientation_map = zeros(size(human_seg));
edge_map = zeros(size(human_seg));
edgeinfo = [];
for k=1:max(max(human_seg))
   segment=human_seg==k;
   [bmap] = seg2bmap(segment==1);
   [xcoord,ycoord]=find(bmap==1);
   [indices]=find(bmap==1);
   dt=bwdist(segment);
   H=fspecial('gaussian',[7 7],1);
   dt_smoothed=imfilter(dt,H);
   [fx,fy] = gradient(dt_smoothed);
   dmap=sqrt(fx.^2+fy.^2);
   fx=fx./dmap;
   fy=fy./dmap;
   tangents=[fy(indices) -fx(indices)];
   tan_angle=atan2(tangents(:,2),tangents(:,1));
   orientation_map(indices)=tan_angle;
   edge_map = edge_map|bmap;
   
%    segment_smoothed=imfilter(double(segment),H);
%    [fx,fy] = gradient(segment_smoothed);
%    dmap=sqrt(fx.^2+fy.^2);
%    [ ~, TO_orientation, temp ] = subpix_TO_correction( dmap,orientation_map,0,bmap);
%    edgeinfo = [edgeinfo;temp];
end
end