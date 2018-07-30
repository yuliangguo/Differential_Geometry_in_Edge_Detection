close all; clear; clc;
addpath(genpath('toolbox'));

% path
inpath = 'data/';
outpath = 'plots/MDCA/';
mkdir(outpath);
files = dir([inpath '*.png']);

Cvrvelets_diff_avg = [];
Cvrvelets_diff_max = [];
MDCA_diff_avg = [];
MDCA_diff_max = [];

for f = 1:length(files)
    name = files(f).name(1:end-4);
    fprintf(1,'Image %d: %s',f,name);

    % load curvelets of each edge
    [ chain, config, edg, opts] = load_cvlet([inpath name '.cvlet']);
    
    % load ground-truth: x_gt, y_gt, k_gt
    load([inpath name '_gt.mat'])
    
    % load MDCA estimate
    load([inpath name '_MDCA.mat'])
    num_points = length(k);
    k_gt_temp = ones(num_points,1);
    for p = 1:length(k_gt_temp),
        dist = (x(p) - x_gt).^2 + (y(p) - y_gt).^2;
        [mindist,gt_ind] = min(dist);
        k_gt_temp(p) = k_gt(gt_ind);
        if(k_gt_temp(p)<0), k(p) = -k(p); end
    end
    
    % match estimate values with GT
    [kdata,kgtdata] = get_curvature_data_points( x, y, k_gt_temp, ...
        edg, chain, config);
    
    % plot the the test results
    figure(1); clf;
    pointid = 1:num_points;
    plot(pointid(~isnan(kgtdata)),kdata(~isnan(kgtdata)),'r', 'LineWidth',2); hold on;
    plot(pointid,k,'b', 'LineWidth',2);
    plot(pointid,k_gt_temp,'g--', 'LineWidth',2);
    ylim([-0.3,0.3]);
    xlim([0,num_points]);
    set(gca, 'FontSize', 22);
    legend('Curvelets','MDCA','Ground Truth','Location','NorthWest');
%     title(['Curvature Estimation of ' ])
    print(gcf,'-dpdf',[outpath 'curvature_estimate_of_' name '.pdf']);
    cmd = ['!pdfcrop ' [outpath 'curvature_estimate_of_' name '.pdf'] ' ' [outpath 'curvature_estimate_of_' name '.pdf'] ];
    eval(cmd);
    
    kdata_smooth = smooth(repmat(kdata(~isnan(kgtdata)),3,1),...
        2*round(num_points/(max(abs(kdata))*1000)));
    kdata_smooth = kdata_smooth(length(kdata_smooth)/3+1:2*length(kdata_smooth)/3);
    figure(2); clf;
    plot(pointid(~isnan(kgtdata)),kdata_smooth,'r', 'LineWidth',2); hold on;
    plot(pointid,k,'b', 'LineWidth',2);
    plot(pointid,k_gt_temp,'g--', 'LineWidth',2); 
    ylim([-0.3,0.3]);
    xlim([0,num_points]);
    set(gca, 'FontSize', 22);
    legend('Curvelets','MDCA','Ground Truth','Location','NorthWest');
%     title(['Curvature Estimation of ' ])
    print(gcf,'-dpdf',[outpath 'curvature_estimate_of_' name '_smooth.pdf']);
    cmd = ['!pdfcrop ' [outpath 'curvature_estimate_of_' name '_smooth.pdf'] ' ' [outpath 'curvature_estimate_of_' name '_smooth.pdf'] ];
    eval(cmd);
    
    Curvelets_diff = abs(kgtdata(~isnan(kgtdata)) -  kdata_smooth);
    MDCA_diff = abs(kgtdata(~isnan(kgtdata)) -  k(~isnan(kgtdata)));
    
    Cvrvelets_diff_avg = [Cvrvelets_diff_avg mean(Curvelets_diff)];
    Cvrvelets_diff_max = [Cvrvelets_diff_max max(Curvelets_diff)];
    MDCA_diff_avg = [MDCA_diff_avg mean(MDCA_diff)];
    MDCA_diff_max = [MDCA_diff_max max(MDCA_diff)];
    
end