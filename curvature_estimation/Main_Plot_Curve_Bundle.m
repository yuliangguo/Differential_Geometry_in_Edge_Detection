col = {[0.75 0.9 1], [1 0.75 0.9], [1 0.9 0.75], [0.9 1 0.75]};
% load('bundle.mat');
% kth = 0.15;
% figure(1);
% cum_bmap = ones(size( Kmin{1}));
% for i = 1:length(Kmin)
%     figure(1);
%     if i~=1
%         hold on;
%         fill(bundle_polygon(:,1),bundle_polygon(:,2),col{col_ind});
%         hold off;
%     end
%     bmap = Kmin{order(i)}<kth & Kmax{order(i)}>kth;
%     dbundle_polygon = regionprops(bmap,'Extrema');
%     dbundle_polygon = dbundle_polygon.Extrema;
%     figure(2), imagesc(bmap);
%     cum_bmap = cum_bmap & bmap;
%     bundle_polygon = regionprops(cum_bmap,'Extrema');
%     bundle_polygon = bundle_polygon.Extrema;
%     figure(1);
%     col_ind = randi(4,1); 
%     hold on;
%     fill(dbundle_polygon(:,1),dbundle_polygon(:,2),col{col_ind});
%     fill(bundle_polygon(:,1),bundle_polygon(:,2),'r');
%     hold off;
%     
%     
%     figure(3);
%     if i~=1
%         hold on;
%         fill(bundle{i-1,2}(:,1),bundle{i-1,2}(:,2),col{col_ind});
%         hold off;
%     end
%     col_ind = randi(4,1); 
%     hold on;
%     fill(bundle{i,1}(:,1),bundle{i,1}(:,2),col{col_ind});
%     fill(bundle{i,2}(:,1),bundle{i,2}(:,2),'r');
%     hold off;
% end
% axis square

close all;
dst_path = 'plots/'
load('data/curve_bundle.mat');
figure(3);
a = gca;
a.FontSize = 18;
xlabel('\delta x (pixels)')
ylabel('\delta \theta (degrees)')
axis([-1.5 1.5 -20 20]);
transformy = @(bundle) [bundle(:,1)*180/pi bundle(:,2)];
bundle(:,1) = cellfun(transformy, bundle(:,1),'UniformOutput', false);
bundle(:,2) = cellfun(transformy, bundle(:,2),'UniformOutput', false);
for i = 1:size(bundle,1)
    figure(3);
    if i~=1
        hold on;
        fill(bundle{i-1,2}(:,2),bundle{i-1,2}(:,1),col{col_ind}, ...
            'EdgeColor','none','LineWidth',0.1);
        hold off;
    end
    col_ind = randi(4,1); 
    hold on;
    fill(bundle{i,1}(:,2),bundle{i,1}(:,1),col{col_ind});
    fill(bundle{i,2}(:,2),bundle{i,2}(:,1),'r','EdgeColor','none',...
        'LineWidth',3);
    hold off;
    print(gcf,'-dpdf',[dst_path 'curve_bundle' num2str(i) '.pdf'])
end
