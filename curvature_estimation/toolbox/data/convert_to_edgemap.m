function [ edgemap, edgeid, thetamap ] = convert_to_edgemap( edgeinfo, h, w )
edgemap = zeros(h,w);
thetamap = zeros(h,w);
edgeid = cell(h,w);
for k = 1:size(edgeinfo,1)
    xpos = min([w,max([1,round(edgeinfo(k,1)+1)])]);
    ypos = min([h,max([1,round(edgeinfo(k,2)+1)])]);
    edgemap(ypos,xpos) = max(edgemap(ypos,xpos),edgeinfo(k,4));
    thetamap(ypos,xpos) = edgeinfo(k,3);
    edgeid{ypos,xpos} = [edgeid{ypos,xpos}, k];
end
end

