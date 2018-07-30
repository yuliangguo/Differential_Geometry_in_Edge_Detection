close all; clear; clc;
addpath(genpath('toolbox'));

% path
inpath = 'data/';
outpath = 'plots/';
mkdir(outpath);

cvletfile = 'ellipse.cvlet';
cemfile = 'ellipse.cem';

% load curvelets of each edge
[ chain, config, edg, opts] = load_cvlet([inpath cvletfile]);

% load the contour of each ellipse to find out the ellipse an edge belongs
% to.
cem = load_contours([inpath cemfile]);

%% ellipse
% ellipse parameters for gt curvature
contour_ind = [1;5;3;4;6;2];
f1 = [20 40; 176 38; 40 131; 50 220; 227 150; 33 310];
f2 = [80 40; 263 67; 177 131; 144 200; 232 353; 152 338];
e = [0.8; 0.6; 0.95; 0.75; 0.85; 0.55];

w = atan2(f2(:,2)-f1(:,2),f2(:,1)-f1(:,1));
c = (f1+f2)/2;
a = 1/2*sqrt((f2(:,1)-f1(:,1)).^2+(f2(:,2)-f1(:,2)).^2);
b = a.*sqrt(1-e(:).^2);


curvenum = length(cem{2});
kdata = cell(curvenum,1);
kmaxdata = cell(curvenum,1);
kmindata = cell(curvenum,1);
k_norm_error_data = cell(curvenum,1);
kgt = cell(curvenum,1);
vertex_k = zeros(curvenum,1);
vertex_k_gt = zeros(curvenum,1);
co_vertex_k = zeros(curvenum,1);
co_vertex_k_gt = zeros(curvenum,1);

% record curveture data
for i = 1:curvenum
    x = @(t) a(i)*cos(t)*cos(w(i)) - b(i)*sin(t)*sin(w(i)) + c(i,1);
    y = @(t) a(i)*cos(t)*sin(w(i)) + b(i)*sin(t)*cos(w(i)) + c(i,2);
    xt = @(t) -a(i)*sin(t)*cos(w(i)) - b(i)*cos(t)*sin(w(i));
    yt = @(t) b(i)*cos(t)*cos(w(i)) - a(i)*sin(t)*sin(w(i));
    xtt = @(t) -a(i)*cos(t)*cos(w(i)) + b(i)*sin(t)*sin(w(i));
    ytt = @(t) -b(i)*sin(t)*cos(w(i)) - a(i)*cos(t)*sin(w(i));
    
    sample_num = ceil(3*pi*a(i));
    [kdata{i},kgt{i},kmaxdata{i},kmindata{i}] = get_curvature_data( xtt, ...
        ytt, xt, yt, x, y, cem{2}{contour_ind(i)}, chain, config,sample_num );
    
    % smooth the data
    num_points = length(kgt{i});
    kdata{i} = smooth(repmat(kdata{i}(~isnan(kgt{i})),3,1),...
        3*round(num_points/(max(abs(kdata{i}))*1000)));
    kdata{i} = kdata{i}(length(kdata{i})/3+1:2*length(kdata{i})/3);

    k_norm_error_data{i} = (kdata{i}-kgt{i})./kgt{i}; 

    num_samples = length(kdata{i});
    [vertex_k_gt(i), ind] = max(kgt{i});
    vertex_k(i) = kdata{i}(ind);
    [co_vertex_k_gt(i), ind] = min(kgt{i});
    co_vertex_k(i) = kdata{i}(ind);
end

[vertex_k_gt, ind] = sort(vertex_k_gt, 'ascend');
vertex_k = vertex_k(ind);
[co_vertex_k_gt, ind] = sort(co_vertex_k_gt, 'ascend');
co_vertex_k =  co_vertex_k(ind);

% plot the the test results
for i = 1:curvenum
    clf
    %% plot estimated vs. GT
    figure(1);
    plot(1:length(kdata{i}),kdata{i},'r.');hold on
    plot(1:length(kgt{i}),kgt{i},'g'); 
%     plot(1:length(kmindata{i}),kmindata{i},'k');
%     plot(1:length(kmaxdata{i}),kmaxdata{i},'b');
    hold off;
    ylim([-0.2,0.2]);
    xlim([0,max(length(kdata{i}), length(kgt{i}))]);
    set(gca, 'FontSize', 18);
%     legend('Estimated Curvature','Ground Truth','Minimum Curvature',...
%         'Maximum Curvature','Location','SouthEest');
    legend('Estimated Curvature','Ground Truth','Location','SouthEast');
    title(['Curvature Estimation of Curve ' num2str(contour_ind(i))])
    print(gcf,'-dpdf',[outpath 'curvature_estimate_of_ellipse_' ...
        num2str(contour_ind(i))]);
    
    %% plot the Normalized Curvature Difference
    figure(2);
    plot(1:length(k_norm_error_data{i}),k_norm_error_data{i},'r'); 
%     plot(1:length(kmindata{i}),kmindata{i},'k');
%     plot(1:length(kmaxdata{i}),kmaxdata{i},'b');hold off;
%     ylim([0,0.5]);
    xlim([0,length(k_norm_error_data{i})]);
    set(gca, 'FontSize', 15);
%     legend('Estimated Curvature','Ground Truth','Minimum Curvature',...
%         'Maximum Curvature','Location','SouthEest');
    legend('Normalized Curvature Difference','Location','NorthEast');
    title(['Normalized Curvature Difference ' num2str(contour_ind(i))])
    print(gcf,'-dpdf',[outpath 'curvature_norm_diff_of_ellipse_' ...
        num2str(contour_ind(i))]);
end

%% scatter of curvature estimation all ellipses
h=figure(3);hold on;

for i = 1:curvenum
    plot(kgt{i}, kdata{i},'r.')
    plot(kgt{i}, kgt{i},'g--')
end
hold off;
set(gca,'FontSize',18)

% axis([0.02 vertex_k_gt(end) 0.02 vertex_k_gt(end)])
grid on;
grid minor;
legend('Estimated Curvature',...
    'Ground Truth', 'Location','SouthEast')
xlabel('Ground Truth Curvature');
ylabel('Estimated Curvature')
title('Curvature Estimates of Ellipses')
print( h, '-dpdf', [outpath 'curvature_estimates_of_ellipses_all']);

% set(gca,'XScale','log');
% set(gca,'YScale','log');
% print( h, '-dpdf', [outpath 'curvature_estimates_of_ellipses_all_loglog']);
%% scatter of normalized curvature diff all ellipses
h=figure(4);hold on;

for i = 1:curvenum
    plot(kgt{i}, abs((kdata{i}-kgt{i})./kgt{i}),'r.')
end
hold off;
set(gca,'FontSize',18)

% axis([0.02 vertex_k_gt(end) 0.02 vertex_k_gt(end)])
grid on;
grid minor;
legend('Normalized Curvature Difference',...
     'Location','NorthWest')
xlabel('Ground Truth Curvature');
ylabel('Normalized Curvature Difference')
title('Normalized Curvature Difference of Ellipses')
print( h, '-dpdf', [outpath 'norm_curvature_diff_ellipses_all']);

% set(gca,'XScale','log');
% set(gca,'YScale','log');
% print( h, '-dpdf', [outpath 'norm_curvature_diff_ellipses_all_loglog']);

%% plot the average estimate values at the Vertex
h=figure(3);
plot(vertex_k_gt, vertex_k,'rs-','LineWidth',3)
hold on;
plot(vertex_k_gt,vertex_k_gt,'gs-','LineWidth',3);
hold off;
set(gca,'FontSize',15)

axis([0.02 vertex_k_gt(end) 0.02 vertex_k_gt(end)])
grid on;
grid minor;
legend('Estimated Curvature',...
    'Ground Truth', 'Location','NorthWest')
xlabel('Ground Truth Curvature');
ylabel('Estimated Curvature')
title('Curvature Estimates at Ellipse Vertex')
print( h, '-dpdf', [outpath 'curvature_estimates_of_ellipses_vertex']);

set(gca,'XScale','log');
set(gca,'YScale','log');
print( h, '-dpdf', [outpath 'curvature_estimates_of_ellipses_vertex_loglog']);

%% plot the norm diff at the Vertex
h=figure(3);
plot(vertex_k_gt, abs(vertex_k-vertex_k_gt)./vertex_k_gt,'rs-','LineWidth',3)


xlim([0.02 vertex_k_gt(end)])
ylim([0.1 0.4])
grid on;
grid minor;
legend('Normalized Curvature Difference',...
    'Location','NorthWest')
xlabel('Ground Truth Curvature');
ylabel('Normalized Curvature Difference')
title('Normalized Curvature Difference at Ellipse Vertex')
set(gca,'FontSize',15)

print( h, '-dpdf', [outpath 'curvature_norm_error_ellipses_vertex']);

set(gca,'XScale','log');
set(gca,'YScale','log');
print( h, '-dpdf', [outpath 'curvature_norm_error_ellipses_vertex_loglog']);

% %% plot the average estimate values at the Co-Vertex
% h=figure(4);
% plot(co_vertex_k_gt, co_vertex_k,'rs-','LineWidth',3)
% hold on;
% plot(co_vertex_k_gt, co_vertex_k_gt,'gs-','LineWidth',3);
% hold off;
% set(gca,'FontSize',15)
% 
% % axis([0.01 GT(end) 0.01 GT(end)])
% grid on;
% grid minor;
% legend('Estimated Curvature',...
%     'Ground Truth', 'Location','NorthWest')
% xlabel('Ground Truth Curvature');
% ylabel('Estimated Curvature')
% title('Curvature Estimates at Ellipse Co-Vertex')
% print( h, '-dpdf', [outpath 'curvature_estimates_of_ellipses_co-vertex']);
% 
% set(gca,'XScale','log');
% set(gca,'YScale','log');
% print( h, '-dpdf', [outpath 'curvature_estimates_of_ellipses_vertex_loglog']);



%% plane curve 1
cvletfile = 'plane_curve1.cvlet';
cemfile = 'plane_curve1.cem';

% load curvelets of each edge
[ chain, config, edg, opts] = load_cvlet([inpath cvletfile]);

% load the contour of each ellipse to find out the ellipse an edge belongs
% to.
cem = load_contours([inpath cemfile]);

dx = 160;
dy = 200;
amp = 50;
contour_ind = [1,2,3,6,7];

curvenum = length(contour_ind);
kdata = cell(curvenum,1);
kmaxdata = cell(curvenum,1);
kmindata = cell(curvenum,1);
kgt = cell(curvenum,1);
k_norm_error_data = cell(curvenum,1);
% record curveture data
for i = 1:curvenum
    % plane curve 1
    x = @(t) amp*(sin(2.^t) - 1.7).*cos(t)+dx;
    y = @(t) amp*(sin(2.^t) - 1.7).*sin(t)+dy;
    xt = @(t) -amp*(sin(2.^t) - 1.7).*sin(t) + amp*cos(2.^t)*2.^t*log(2).*cos(t);
    yt = @(t) amp*(sin(2.^t) - 1.7).*cos(t)  + amp*cos(2.^t)*2.^t*log(2).*sin(t);
    xtt = @(t) -amp*(sin(2.^t) - 1.7).*cos(t) - amp*cos(2.^t)*2.^t*log(2).*sin(t)- ...
        amp*cos(2.^t)*2.^t*log(2).*sin(t) - amp*sin(2.^t)*4.^t*2*log(2).*cos(t) + ...
        amp*cos(2.^t)*2.^t*2*log(2).*cos(t);
    ytt = @(t) -amp*(sin(2.^t) - 1.7).*sin(t) + amp*cos(2.^t)*2.^t*log(2).*cos(t)+ ...
        amp*cos(2.^t)*2.^t*log(2).*cos(t) - amp*sin(2.^t)*4.^t*2*log(2).*sin(t) + ...
        amp*cos(2.^t)*2.^t*2*log(2).*sin(t);
    sample_num = 600;
    [kdata{i},kgt{i},kmaxdata{i},kmindata{i}] = get_curvature_data( xtt, ...
        ytt, xt, yt, x, y, cem{2}{contour_ind(i)}, chain, ...
        config,sample_num );
    
    % smooth the data
    num_points = length(kgt{i});
    kdata{i} = smooth(repmat(kdata{i}(~isnan(kgt{i})),3,1),...
        3*ceil(num_points/(max(abs(kdata{i}))*1000)));
    kdata{i} = kdata{i}(length(kdata{i})/3+1:2*length(kdata{i})/3); 
    
    k_norm_error_data{i} = abs(kdata{i}-kgt{i})./abs(kgt{i}); 

end

% plot the the test results
for i = 1:curvenum
    clf
    if(isempty(kdata{i}))
        continue;
    end
    plot(1:length(kdata{i}),kdata{i},'r'); hold on;
    plot(1:length(kgt{i}),kgt{i},'g'); 
%     plot(1:length(kmindata{i}),kmindata{i},'k');
%     plot(1:length(kmaxdata{i}),kmaxdata{i},'b');
    hold off;
    ylim([-0.3,0.3]);
    set(gca, 'FontSize', 18);
    legend('Estimated Curvature','Ground Truth','Minimum Curvature',...
        'Maximum Curvature','Location','NorthWest');
%     title(['Curvature Estimation of Curve ' num2str(contour_ind(i))])
    print(gcf,'-dpdf',[outpath 'curvature_estimate_of_plane_curve_1_' ...
        num2str(contour_ind(i))]);
    
    %% plot the Normalized Curvature Difference
    figure(2);
    plot(1:length(k_norm_error_data{i}),k_norm_error_data{i},'r'); 
%     plot(1:length(kmindata{i}),kmindata{i},'k');
%     plot(1:length(kmaxdata{i}),kmaxdata{i},'b');hold off;
%     ylim([0,0.5]);
    xlim([0,length(k_norm_error_data{i})]);
    set(gca, 'FontSize', 15);
%     legend('Estimated Curvature','Ground Truth','Minimum Curvature',...
%         'Maximum Curvature','Location','SouthEest');
    legend('Normalized Curvature Difference','Location','NorthEast');
%     title(['Normalized Curvature Difference ' num2str(contour_ind(i))])
    print(gcf,'-dpdf',[outpath 'curvature_norm_diff_of_plane_curve_1_' ...
        num2str(contour_ind(i))]);
end

%% scatter of curvature estimation all samples in the plane shape
close all;
h=figure(3);hold on;

for i = 2:curvenum
    plot(kgt{i}, kdata{i},'r.')
    plot(kgt{i}, kgt{i},'g--')
end
hold off;
axis([-0.2 0.2 -0.2 0.2])

% axis([0.02 vertex_k_gt(end) 0.02 vertex_k_gt(end)])
grid on;
grid minor;
legend('Estimated Curvature',...
    'Ground Truth', 'Location','SouthEast')
xlabel('Ground Truth Curvature');
ylabel('Estimated Curvature')
title('Curvature Estimates of Random Shape')
set(gca,'FontSize',15)
print( h, '-dpdf', [outpath 'curvature_estimates_of_Random_Shape']);

% set(gca,'XScale','log');
% set(gca,'YScale','log');
% print( h, '-dpdf', [outpath 'curvature_estimates_of_ellipses_all_loglog']);
%% scatter of normalized curvature diff all samples in the plane shape
% clear all;
h=figure(4);hold on;

for i = 2:curvenum
    plot(kgt{i}, abs((kdata{i}-kgt{i})./kgt{i}),'r.')
end
hold off;

grid on;
grid minor;
legend('Normalized Curvature Difference',...
     'Location','NorthWest')
xlabel('Ground Truth Curvature');
ylabel('Normalized Curvature Difference')
title('Normalized Curvature Difference of Random Shape')
axis([-0.2 0.2 0 100])
set(gca,'FontSize',14)

print( h, '-dpdf', [outpath 'norm_curvature_diff_Random_Shape']);

% set(gca,'XScale','log');
% set(gca,'YScale','log');
% print( h, '-dpdf', [outpath 'norm_curvature_diff_ellipses_all_loglog']);

