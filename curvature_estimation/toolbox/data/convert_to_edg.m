function [ edginfo ] = convert_to_edg( edgemap, thetamap )
[r,c] = find(edgemap~=0);
ind = sub2ind(size(edgemap),r,c);
edginfo = [c-1,r-1,thetamap(ind),edgemap(ind)];
end

