function [ curvefrag, cirarcfrag ] = convert_cfrag( chain, edgeinfo, curveinfo )
% Input:
%   chain: edge id chains of curves
%   curveinfo: [isForward ref_pt.x() ref_pt.y() ref_theta pt.x() pt.y() 
%       theta k length property center.x() center.y() ]
%   edgeinfo: properties of each edges
% Output:
%   curvefrag: {1 x numcurve} each cell stores the edgeinfo belongs to that
%   curve
%   cirarcfrag: [numcurve x 8] each row stores the underlying circular arc
%   position and properties

numcurve = size(chain,1);
curvefrag = cell(1,numcurve);
cirarc_output = false;
if nargin>2 && nargout>1
    if size(curveinfo,2)==12
        cirarcfrag = zeros(numcurve,8);
        cirarc_output = true;
    else
        warning('Only one output is expected.');
        cirarcfrag = [];
    end
end
for i = 1:numcurve
    ind = chain(i,:);
    ind = ind(ind~=0);
    numedge = numel(ind);
    curvefrag{i} = zeros(numedge,5);
    curvefrag{i} = edgeinfo(ind,:);
    
    if cirarc_output
        % start point
        cirarcfrag(i,1:3) = curveinfo(i,5:7);
        % center
        cirarcfrag(i,4:5) = curveinfo(i,11:12);
        % curvature
        cirarcfrag(i,6) = curveinfo(i,8);
        % length
        cirarcfrag(i,7) = curveinfo(i,9);
        % quality
        cirarcfrag(i,8) = curveinfo(i,10);
    end
end

end

