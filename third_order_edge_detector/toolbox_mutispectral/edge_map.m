function edge_map( X, Y, ANG, length, color )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% Note: X:1*n
%       Y:1*n
hold on;
for p=1:size(X,2)
    plot([(X(p)+length*cos(ANG(p))) (X(p)-length*cos(ANG(p)))]...
        ,[(Y(p)+length*sin(ANG(p))) (Y(p)-length*sin(ANG(p)))], ...
        color);
end
hold off;

end

