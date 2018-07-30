function edges_new = update_edges_via_curvelets(edges, chain,info)

% update each edge's location with the best estimated from curvelets
% forms curvelets
edges_new = edges;

for i = 1:size(edges, 1)
    idx = find(chain(:,1) == i);
%     anchor = mean(info(idx, 1:2));
    if(length(idx==1))
        edge = mean(info(idx, 5:7));
        edges_new(i,1:3) = edge;
    end
end

% edges_new(:,1:2) = edges_new(:,1:2)-1;