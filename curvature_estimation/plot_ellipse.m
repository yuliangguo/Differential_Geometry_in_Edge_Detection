function points = plot_ellipse( f1, f2, e, col)
% Input:
%       f1, f2: [1x2] focal points;
%       e: the eccentricity of an ellipse, f/a;
if nargin<4
    col = 'r';
end
a = 1/2*sqrt((f2(1)-f1(1))^2+(f2(2)-f1(2))^2);
b = a*sqrt(1-e^2);
% approximate circumference
c = pi*(3*(a+b)-sqrt((3*a+b)*(a+3*b)));
% disp(c)
t = linspace(0,2*pi,ceil(c));
X = a*cos(t);
Y = b*sin(t);
w = atan2(f2(2)-f1(2),f2(1)-f1(1));
x = (f1(1)+f2(1))/2 + X*cos(w) - Y*sin(w);
y = (f1(2)+f2(2))/2 + X*sin(w) + Y*cos(w);
hold on;
plot(x,y,[col,'.']);
hold off;
points = [x;y;X;Y;t];
end

