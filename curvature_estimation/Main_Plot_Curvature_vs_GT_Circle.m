close all; clear; clc;
addpath(genpath('toolbox'));

% path
inpath = 'data/';
outpath = 'plots/';
mkdir(outpath);
cvletfile = 'circle_fill_blur.cvlet';

cemfile = 'circle_fill_blur.cem';

% load curvelets of each edge
[ chain, config, edg, opts] = load_cvlet([inpath cvletfile]);
% load the contour of each circle to find out the curvature ground truth 
% corresponding to the edge.
cem = load_contours([inpath cemfile]);

curvenum = length(cem{2});
curvatureavg = zeros(1,curvenum);
curvaturemax = zeros(1,curvenum);
curvaturemin = zeros(1,curvenum);
% radius = zeros(1,curvenum);
GT = [1/75, 1/35, 1/10, 1/15, 1/1.5, 1/50, 1/2.5, 1/20, 1/100, ...
    1/3.5, 1/25, 1/5, 0];
kdata = cell(curvenum,1);
kmaxdata = cell(curvenum,1);
kmindata = cell(curvenum,1);
k_norm_error_data = cell(curvenum,1);
cvlet_min_len = 2;
% record curveture data
for i = 1:curvenum
    count = 0;
    num = 0;
    cumk = 0;
    cvmax = 0;
    cvmin = 100;
    cvkdata = [];
    cvkmaxdata = [];
    cvkmindata = [];
    for j = 1:size(cem{2}{i},1)
        edgeID = cem{2}{i}(j,6);
        cvlets_ind = find(chain(:,1)==edgeID);
        cvlets_chain = chain(cvlets_ind,:);
        cvlets_config = config(cvlets_ind,:);
        for k = 1:length(cvlets_ind) 
            % iterate all the curvelets ankered at this point
            num  = num+1;
            if all(cvlets_chain(k,:)) || find(cvlets_chain(k,:)==0,1)>cvlet_min_len+1
                count = count+1;
                tempk = abs(cvlets_config(k,8));
                cvkdata = [cvkdata cvlets_config(k,8)];
                cvkmaxdata = [cvkmaxdata cvlets_config(k,2)];
                cvkmindata = [cvkmindata cvlets_config(k,3)];
                cumk = cumk + tempk;
%                 disp(1/tempk);
%                 if tempk>3*10^-3
%                     cumr = cumr + 1/tempk;
%                     count = count + 1;
%                 end
                
                if(tempk>cvmax)
                    cvmax = tempk;
                end
                if(tempk<cvmin)
                    cvmin = tempk;
                end
            end
        end
    end
    kdata{i} = cvkdata;
    
    % smooth the data
    num_points = length(kdata{i});
%     kdata{i} = smooth(repmat(kdata{i},3,1),...
%         3*ceil(num_points/(max(abs(kdata{i}))*1000)));
    kdata{i} = smooth(repmat(kdata{i},3,1),...
        10);
    kdata{i} = kdata{i}(length(kdata{i})/3+1:2*length(kdata{i})/3);
    
    k_norm_error_data{i} = abs(cvkdata-GT(i))/GT(i); 
    kmaxdata{i} = cvkmaxdata;
    kmindata{i} = cvkmindata;
    curvatureavg(i) = cumk/num;
%     radius(i) = cumr/count;
    disp(['data:' num2str(num) ', used:' num2str(count)...
        ', percent:' num2str(count/num)])
    curvaturemax(i) = cvmax;
    curvaturemin(i) = cvmin;
end

[GT,ind] = sort(GT,'ascend');
% GT = [0, 0.01, 0.0133, 0.02, 0.0286, 0.04, 0.05, 0.0667, 0.1, 0.2,
% 0.2857, 0.4, 0.6667]
% the last three values exceed the upper bound of estimation value
GT = GT(1:end-3);
curvatureavg = curvatureavg(ind(1:end-3));
curvaturemax = curvaturemax(ind(1:end-3));
curvaturemin = curvaturemin(ind(1:end-3));
% for radius, remove the infinite gt value
radius = 1./curvatureavg(2:end);
radiusGT = 1./GT(2:end);

kdata = kdata(ind(1:end-3));
kmaxdata = kmaxdata(ind(1:end-3));
kmindata = kmindata(ind(1:end-3));
k_norm_error_data = k_norm_error_data(ind(1:end-3));


figcount = 1;
% %% histogram of all estimated values
% for i = 1:length(kdata)
%     h = figure(figcount);
%     figcount = figcount+1;
%     [histy, histx] = hist(kdata{i});
%     [histkmaxy, histkmaxx] = hist(kmaxdata{i});
%     [histkminy, histkminx] = hist(kmindata{i});
%     h1 = plot(histx,histy,'r','LineWidth',3); hold on;
%     h2 = plot(histkmaxx,histkmaxy,'b','LineWidth',3); 
%     h3 = plot(histkminx,histkminy,'k','LineWidth',3); 
%     h4 = plot(GT(i)*ones(1,10), linspace(0,1.1*max(histy),10),'g','LineWidth',3); hold off;
%     set(gca,'FontSize',18)
%     xlabel('Estimated Curvature');
%     ylabel('Number of Edges')
%     hl = legend([h1,h2,h3,h4],'Curvature Outputs','$k_{max}$','$k_{min}$','Ground-Truth');
%     set(hl,'Interpreter','latex');
%     print( h, '-dpdf', [outpath 'curvature_output_histogram_id_' num2str(i)]);
% end
% 
% %% Confidence Interval
% 
% for i = 1:length(kdata)
%     h = figure(figcount);
%     clf;
%     figcount = figcount+1;
%     leftinterval = kdata{i} - kmindata{i};
%     rightinterval = kmaxdata{i} - kdata{i};
%     interval = kmaxdata{i} - kmindata{i};
%     [histlefty, histleftx] = hist(leftinterval);
%     [histrighty, histrightx] = hist(rightinterval);
%     [histally, histallx] = hist(interval);
%     h1 = plot(histleftx,histlefty,'r','LineWidth',3); hold on;
%     h2 = plot(histrightx,histrighty,'b','LineWidth',3); 
%     h3 = plot(histallx,histally,'g','LineWidth',3); 
%     set(gca,'FontSize',18)
%     xlabel('Confifence Interval of Estimated Curvature');
%     ylabel('Number of Edges')
%     hl = legend([h1,h2,h3],'$\overline{k}-k_{min}$','$k_{max}-\overline{k}$',...
%         '$k_{max}-k_{min}$');
%     set(hl,'Interpreter','latex');
%     print( h, '-dpdf', [outpath 'curvature_data_confidence_interval_id_' num2str(i)]);
% end



%% single data points visualization
% max_num_data = 100;
% 
% for i = 1:length(kdata)
%     h = figure(figcount);
%     clf;
%     figcount = figcount+1;
%     num_data = min([numel(kdata{i}) max_num_data]);
%     ind = randperm(numel(kdata{i}), num_data);
%     hold on;
%     for r = 1:length(ind)
%         h1 = plot([kmindata{i}(ind(r)) kdata{i}(ind(r)) kmaxdata{i}(ind(r))],...
%             ind(r)*ones(1,3),'rs-','LineWidth',1);
%     end
%     h4 = plot(GT(i)*ones(1,10), linspace(0,1.1*numel(kdata{i}),10),'g','LineWidth',3); hold off;
%     set(gca,'FontSize',18)
%     axis([-0.3,0.3,0,numel(kdata{i})])
%     xlabel('Estimated Curvature');
%     ylabel('Curvelets ID')
%     hl = legend([h1,h4],'Curvature Output','Ground-Truth');
% %     set(hl,'Interpreter','latex');
%     print( h, '-dpdf', [outpath 'curvature_data_output_visualization_id_' num2str(i)]);
% end


for i = 1:length(kdata)

    if(isempty(kdata{i}))
       continue; 
    end
    h=figure(1);
    plot(1:length(kdata{i}),kdata{i},'r-');hold on
    plot(1:length(kdata{i}),GT(i)*ones(size(kdata{i})),'g-'); 
%     plot(1:length(kmindata{i}),kmindata{i},'k');
%     plot(1:length(kmaxdata{i}),kmaxdata{i},'b');
    hold off;
    ylim([-0.05,0.35]);
    xlim([0,length(kdata{i})]);
    xlabel('Edge ID')
    ylabel('Estimated Curvature');
    legend('Estimated Curvature','Ground Truth','Location','NorthEast');
    set(gca, 'FontSize', 18);
    
%     legend('Estimated Curvature','Ground Truth','Minimum Curvature',...
%         'Maximum Curvature','Location','SouthEest');
    if(i==1)
        title(['Curvature Estimation of Line'])
    else
        title(['Curvature Estimation of Circle: R=' num2str(radiusGT(i-1))])
    end
    print(h,'-dpdf',[outpath 'curvature_estimate_of_circles_id_' num2str(i)]);
     
end

%% plot the average estimate values with StD
h = figure(figcount); figcount = figcount+1;
meank = cellfun(@mean, kdata);
stdk = cellfun(@std, kdata);
h1 = errorbar(GT',meank,stdk,'r','LineWidth',3);
hold on;
h4 = plot(GT,GT,'g--','LineWidth',3);
hold off;
set(gca,'FontSize',15)

axis equal;
axis([0.01 GT(end) 0.01 meank(end)+stdk(end)])
grid on;
grid minor;
legend([h1 h4], 'Estimated Curvature',...
    'Ground Truth', 'Location','NorthWest')
xlabel('Ground Truth Curvature');
ylabel('Estimated Curvature')

print( h, '-dpdf', [outpath 'curvature_errorbar_output']);

set(gca,'XScale','log');
set(gca,'YScale','log');
print( h, '-dpdf', [outpath 'curvature_errorbar_output_loglog']);

%% plot the Normalized Curvature Difference vs. GT curvature
h = figure(figcount); figcount = figcount+1;
meank = cellfun(@mean, k_norm_error_data);
% stdk = cellfun(@std, k_dist_norm_data);
% h1 = errorbar(GT(2:end)',meank(2:end),stdk(2:end),'r','LineWidth',3);
h1 = plot(GT(2:end), meank(2:end),'r-','LineWidth',3);

% hold on;
% h4 = plot(GT,GT,'g--','LineWidth',3);
% hold off;
set(gca,'FontSize',16)
% axis equal;
xlim([0.01 GT(end)])
% axis([0.01 GT(end) 0.01 GT(end)])
grid on;
grid minor;
legend([h1 h4], 'Normalized Curvature Difference',...
     'Location','NorthEast')
xlabel('Ground Truth Curvature');
ylabel('Normalized Curvature Difference')

print( h, '-dpdf', [outpath 'curvature_norm_error_output']);

set(gca,'XScale','log');
set(gca,'YScale','log');
print( h, '-dpdf', [outpath 'curvature_norm_error_output_loglog']);

% %% plot average curvature and GT
% h = figure(figcount); figcount = figcount+1;
% h1 = plot(GT, curvatureavg,'ro-','LineWidth',3);
% hold on;
% h3 = plot(GT,curvaturemax,'k');
% h2 = plot(GT,curvaturemin,'b');
% h4 = plot(GT,GT,'g--','LineWidth',3);
% hold off;
% axis equal;
% axis([GT(1) GT(end) GT(1) GT(end)])
% grid on;
% grid minor;
% set(gca,'FontSize',18)
% legend([h1 h2 h3 h4], 'Estimated Curvature', 'Minimum Curvature', 'Maximum Curvature',...
%     'Ground Truth', 'Location','NorthWest')
% xlabel('Ground Truth Curvature');
% ylabel('Estimated Curvature')
% print( h, '-dpdf', [outpath 'curvature_output']);
% 
% % curvature in loglog space
% h = figure(figcount); figcount = figcount+1;
% h1 = loglog(GT, curvatureavg,'ro-','LineWidth',3);
% hold on;
% % h2 = loglog(GT,curvaturemax),'k');
% % h3 = loglog(GT,curvaturemin),'b');
% h4 = loglog(GT,GT,'g--','LineWidth',3);
% hold off;
% axis equal;
% axis([GT(1) GT(end) GT(1) GT(end)])
% grid on;
% grid minor;
% set(gca,'FontSize',18)
% legend([h1 h4], 'Estimated Curvature','Ground Truth', 'Location','NorthWest')
% xlabel('Ground Truth Curvature');
% ylabel('Estimated Curvature')
% print( h, '-dpdf',[outpath 'curvature_output_loglog']);
% 
% % curvature in semilogx space
% h = figure(figcount); figcount = figcount+1;
% h1 = semilogx(GT, curvatureavg,'ro-','LineWidth',3);
% hold on;
% % h2 = semilogx(GT,curvaturemax,'k');
% % h3 = semilogx(GT,curvaturemin,'b');
% h4 = semilogx(GT,GT,'g--','LineWidth',3);
% hold off;
% axis equal;
% axis([GT(1) GT(end) GT(1) GT(end)])
% grid on;
% grid minor;
% set(gca,'FontSize',18)
% legend([h1 h4], 'Estimated Curvature','Ground Truth', 'Location','NorthWest')
% xlabel('Ground Truth Curvature');
% ylabel('Estimated Curvature')
% print( h, '-dpdf',[outpath 'curvature_output_semilogx']);
% 
% % curvature in semilogy space
% h = figure(figcount); figcount = figcount+1;
% h1 = semilogy(GT, curvatureavg,'ro-','LineWidth',3);
% hold on;
% % h2 = semilogy(GT, curvaturemax, 'k');
% % h3 = semilogy(GT, curvaturemin, 'b');
% h4 = semilogy(GT, GT, 'g--', 'LineWidth', 3);
% hold off;
% axis equal;
% axis([GT(1) GT(end) GT(1) GT(end)])
% grid on;
% grid minor;
% set(gca,'FontSize',18)
% legend([h1 h4], 'Estimated Curvature','Ground Truth', 'Location','NorthWest')
% xlabel('Ground Truth Curvature');
% ylabel('Estimated Curvature')
% print( h, '-dpdf', [outpath 'curvature_output_semilogy']);
% 
% %% plot average radius and GT
% h = figure(figcount); figcount = figcount+1;
% h1 = plot(radiusGT,radius,'ro-','LineWidth',3);
% hold on;
% h4 = plot(radiusGT,radiusGT,'g--','LineWidth',3);
% hold off;
% axis equal;
% axis([radiusGT(end) radiusGT(1) radiusGT(end) radiusGT(1)])
% grid on;
% grid minor;
% set(gca,'FontSize',18)
% legend([h1;h4], 'Estimated Radius', ...
%     'Ground Truth', 'Location','NorthWest')
% xlabel('Ground Truth Radius');
% ylabel('Estimated Radius')
% print( h, '-dpdf', [outpath 'radius_output']);
% 
% % radius in loglog space
% h = figure(figcount); figcount = figcount+1;
% h1 = loglog(radiusGT, radius,'ro-','LineWidth',3);
% hold on;
% h4 = loglog(radiusGT,radiusGT,'g--','LineWidth',3);
% hold off;
% axis equal;
% axis([radiusGT(end) radiusGT(1) radiusGT(end) radiusGT(1)])
% grid minor;
% set(gca,'FontSize',18)
% legend([h1;h4], 'Estimated Radius', ...
%     'Ground Truth', 'Location','NorthWest')
% xlabel('Ground Truth Radius');
% ylabel('Estimated Radius')
% print( h, '-dpdf', [outpath 'radius_output_loglog']);
% 
% % radius in semilogx space
% h = figure(figcount); figcount = figcount+1;
% h1 = semilogx(radiusGT, radius,'ro-','LineWidth',3);
% hold on;
% h4 = semilogx(radiusGT,radiusGT,'g--','LineWidth',3);
% hold off;
% axis equal;
% axis([radiusGT(end) radiusGT(1) radiusGT(end) radiusGT(1)])
% grid on;
% grid minor;
% set(gca,'FontSize',18)
% legend([h1;h4], 'Estimated Radius', ...
%     'Ground Truth', 'Location','NorthWest')
% xlabel('Ground Truth Radius');
% ylabel('Estimated Radius')
% print( h, '-dpdf', [outpath 'radius_output_semilogx']);
% 
% % radius in semilogy space
% h = figure(figcount); figcount = figcount+1;
% h1 = semilogy(radiusGT, radius,'ro-','LineWidth',3);
% hold on;
% h4 = semilogy(radiusGT,radiusGT,'g--','LineWidth',3);
% hold off;
% axis equal;
% axis([radiusGT(end) radiusGT(1) radiusGT(end) radiusGT(1)])
% grid on;
% grid minor;
% set(gca,'FontSize',18)
% legend([h1;h4], 'Estimated Radius', ...
%     'Ground Truth', 'Location','NorthWest')
% xlabel('Ground Truth Radius');
% ylabel('Estimated Radius')
% print( h, '-dpdf', [outpath 'radius_output_semilogy']);
