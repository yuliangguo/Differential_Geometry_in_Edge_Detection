close all; clear; clc;
addpath(genpath('../../toolbox'));

% path
inpath = 'data/';
outpath = 'plots/';
mkdir(outpath);
cvletfile = '009129.cvlet';

% load curvelets of each edge
[ chain, config, edg, opts] = load_cvlet([inpath cvletfile]);

kdata = config(:,8);
kmaxdata = config(:,2);
kmindata = config(:,3);

% plot Confidence Interval
h = figure(1);
clf;
leftinterval = kdata - kmindata;
rightinterval = kmaxdata - kdata;
interval = kmaxdata - kmindata;
[histlefty, histleftx] = hist(leftinterval);
[histrighty, histrightx] = hist(rightinterval);
[histally, histallx] = hist(interval);
% h1 = plot(histleftx,histlefty,'r','LineWidth',3); hold on;
h3 = plot(histrightx,histrighty,'m','LineWidth',3); 
% h3 = plot(histallx,histally,'m','LineWidth',3); 
set(gca,'FontSize',22)
xlabel('Confidence Interval of Curvature');
ylabel('Number of Curvelets')
hl = legend(h3, ...%[h1,h2,h3],'$\overline{k}-k_{min}$','$k_{max}-\overline{k}$',...
    '$k_{max}-k_{min}$');
set(hl,'Interpreter','latex');
print( h, '-dpdf', [outpath 'curvature_data_confidence_interval_' cvletfile(1:end-6)]);

