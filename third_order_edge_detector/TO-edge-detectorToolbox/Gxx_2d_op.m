function Gxx_2d = Gxx_2d_op(sigma, dx, dy)
%function Gxx_2d = Gxx_2d_op(sigma, dx, dy)
%
% d^2G/dx^2 operator with operator shift
%
% (c) LEMS, Brown University
% Amir Tamrakar (amir_tamrakar@brown.edu)
% October 2007

% the size of the operator
w = ceil(sigma*4 + max(abs(dx), abs(dy)));
[x,y] = meshgrid(-w:w,-w:w);

d2G_x = ((x-dx).^2-sigma^2).*exp(-(x-dx).^2/(2*sigma^2))/sqrt(2*pi)/sigma^5;
G_y = exp(-(y-dy).^2/(2*sigma^2))/sqrt(2*pi)/sigma;
Gxx_2d = d2G_x.*G_y;

% % for debug only
% figure;
% axis([-w w -w w]);
% imagesc([-w:w],[-w:w], Gxx_2d); colormap gray;
% hold on;
% th=linspace(0,2*pi,64+1);
% pp=[2.5*Ns*cos(th'), 2.5*Ts*sin(th')];
% pp = pp*[cos(theta+pi/2) sin(theta+pi/2); -sin(theta+pi/2) cos(theta+pi/2)];
% %shift 
% pp(:,1) = pp(:,1)+dx;
% pp(:,2) = pp(:,2)+dy;
% plot(pp(:,1), pp(:,2), 'r', [pp(17,1), pp(49,1)], [pp(17,2), pp(49,2)], 'r', [pp(1,1), pp(33,1)], [pp(1,2), pp(33,2)], 'r');
% axis equal;
% hold off;




