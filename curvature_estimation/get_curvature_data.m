function [kdata,kgt,kmaxdata,kmindata] = get_curvature_data( xtt, ytt, xt, yt, x, y, contour,chain, config,sample_num )
%%%%%%%%% output sample order is the same as the order of computed contour edges
%%%%%%%%% it might be better to follow the ground-truth instead? 
count = 0;
kdata = zeros(size(contour,1),1);
kmaxdata = zeros(size(contour,1),1);
kmindata = zeros(size(contour,1),1);
kgt = zeros(size(contour,1),1);
t = linspace(0,2*pi,sample_num);
sx = x(t);
sy = y(t);
cvlet_min_len = 2;
for j = 1:size(contour,1)
    edgeID = contour(j,6);
    ex = contour(j,1);
    ey = contour(j,2);
    %%%% Find the closed gt points to current contour samples
    dist = (sx-ex).^2 + (sy-ey).^2;
    [mindist,ind] = min(dist);
    if mindist>3, continue; end
%     disp(mindist)
    gt = (ytt(t(ind)).*xt(t(ind))-yt(t(ind)).*xtt(t(ind)))./(xt(t(ind)).^2 ...
        + yt(t(ind)).^2).^(3/2);
    cvlets_ind = find(chain(:,1)==edgeID);
    cvlets_chain = chain(cvlets_ind,:);
    cvlets_config = config(cvlets_ind,:);
    cumk = 0;
    cumkmin = 0;
    cumkmax = 0;
    num_cvlet = 0;
    count = count+1;
    for k = 1:length(cvlets_ind) 
        % iterate all the curvelets ankered at this point
        if all(cvlets_chain(k,:)) || find(cvlets_chain(k,:)==0,1)>cvlet_min_len+1 
            cumk = cumk + cvlets_config(k,8);
            cumkmin = cumkmin + cvlets_config(k,3);
            cumkmax = cumkmax + cvlets_config(k,2);
            num_cvlet = num_cvlet+1;
%                 kdata = [kdata cvlets_config(k,8)];
%                 kmaxdata = [kmaxdata cvlets_config(k,2)];
%                 kmindata = [kmindata cvlets_config(k,3)];
%                 kgt = [kgt gt];
        end
    end
    kdata(count) = cumk/num_cvlet;
    kmaxdata(count) = cumkmax/num_cvlet;
    kmindata(count) = cumkmin/num_cvlet;
    kgt(count) = gt;
end
kdata = kdata(1:count);
kmaxdata = kmaxdata(1:count);
kmindata = kmindata(1:count);
kgt = kgt(1:count);
disp(['data:' num2str(size(contour,1)) ', used:' num2str(count)...
    ', percent:' num2str(count/size(contour,1))])

end

