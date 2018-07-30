function [kdata,kgtdata,kmaxdata,kmindata] = get_curvature_data_points( x_gt, y_gt, k_gt, edg,chain, config )
count = 0;
num_points = size(x_gt,1);
kdata = NaN*ones(num_points,1);
kmaxdata = NaN*ones(num_points,1);
kmindata = NaN*ones(num_points,1);
kgtdata = NaN*ones(num_points,1);
cvlet_min_len = 2;
for j = 1:num_points,
    dist = (edg(:,1) - x_gt(j)).^2 + (edg(:,2) - y_gt(j)).^2;
    [mindist,ind] = min(dist);
    if mindist>3, continue; end
%     disp(mindist)
    cvlets_ind = find(chain(:,1)==ind);
    cvlets_chain = chain(cvlets_ind,:);
    cvlets_config = config(cvlets_ind,:);
    cumk = 0;
    cumkmin = 0;
    cumkmax = 0;
    num_cvlet = 0;
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
%                 kgtdata = [kgtdata gt];
        end
    end
    kdata(j) = cumk/num_cvlet;
    kmaxdata(j) = cumkmax/num_cvlet;
    kmindata(j) = cumkmin/num_cvlet;
    kgtdata(j) = k_gt(j);
end

% disp(['data:' num2str(size(contour,1)) ', used:' num2str(count)...
%     ', percent:' num2str(count/size(contour,1))])

end