function [binary_mask, orient_mask]=convert_contours(contours, interpR)

if(nargin<2)
    interpR = 1;
end

coords=contours{1};
orient_mask=zeros(coords(2),coords(1));
binary_mask=zeros(coords(2),coords(1));
cons=contours{2};

for c=1:length(cons)
    points=cons{c};
    if(interpR>1 && size(points,1) >=9)
        points_new = repmat(points, [interpR, 1]);
        points_new(:,1) = interp(points(:,1), interpR);
        points_new(:,2) = interp(points(:,2), interpR);
        points_new(:,3) = interp(wrapToPi(points(:,3)), interpR);
        points = points_new;
    end
    for d=1:size(points,1)
       xcoord=max(round(points(d,2)+1),1);
       ycoord=max(round(points(d,1)+1),1);
       
       if ( xcoord > coords(2) )
           xcoord=coords(2);
       end
       
       if ( ycoord > coords(1) )
           ycoord=coords(1);
       end
       
       orient_mask(xcoord,ycoord)=points(d,3); 
       binary_mask(xcoord,ycoord)=1;       
    end
    
end

binary_mask = (binary_mask>0);