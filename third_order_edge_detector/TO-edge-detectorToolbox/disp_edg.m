function disp_edg(edg, col)
%function disp_edg(edg, col)
%
% Display the edgemap on the current figure in the given color
%   The edgels are displayed as little line segments.
%
% (c) LEMS, Brown University
% Amir Tamrakar (amir_tamrakar@brown.edu)
% October 2007

hold on;

[m,n] = size(edg);
LinAlpha=0.25;           
ColorValue=[1,0,0]; 
% plot the edgels
length=0.5;
for i=1:m,
    plot([edg(i,1)+1+length*cos(edg(i,3)) edg(i,1)+1-length*cos(edg(i,3))], ...
        [edg(i,2)+1+length*sin(edg(i,3)) edg(i,2)+1-length*sin(edg(i,3))], ...
        'Color', col);%'Color',(1.0-LinAlpha)*(1.0-ColorValue),'EraseMode','xor');
    %plot(edg(i,1), edg(i,2), 'r.');
end

hold off;