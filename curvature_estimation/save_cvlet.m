function save_cvlet(filename, chain, cvletinfo, edges, opts)
% Input:
%   filename: the output file name
%   chain: edge id chains of curves
%   cvinfo: [isForward ref_pt.x() ref_pt.y() ref_theta pt.x() pt.y() 
%       theta k length property]
%   edges: properties of each edges
%   opts: a strucure including the parameters of the curvelet
%       opts.w: image width
%       opts.h: image height
%       opts.nrad: neighborhood radius;
%       opts.dx: maximum distance pertubation;
%       opts.dt: maximum orientation pertubation;
%       opts.token_len: token length;
%
% Code written by Xiaoyan Li
% ref_pt.x() ref_pt.y() are enforced to save K_min and K_max !!!!

% write cvlet header
fprintf(1,'saving into .cvlet file\n');
fout = fopen(filename, 'w');
fprintf(fout, 'CVLET_MAP v2.0\n\n');
fprintf(fout, '[BEGIN HEADER]\n');
fprintf(fout, ['WIDTH=', num2str(opts.w), '\n']);
fprintf(fout, ['HEIGHT=', num2str(opts.h), '\n']);
fprintf(fout, ['EDGE_COUNT=' num2str(size(edges,1)), '\n']);
fprintf(fout, 'CM TYPE=CC2\n');
fprintf(fout, ['N_RADIUS=' num2str(opts.nrad) '\n']);
fprintf(fout, ['DX=' num2str(opts.dx) '\n']);
fprintf(fout, ['DT=' num2str(opts.dt) '\n']);
fprintf(fout, ['TOKEN_LENGTH=' num2str(opts.token_len) '\n']);
fprintf(fout, ['MAX_SIZE_TO_GROUP=' num2str(opts.max_size_to_group) '\n']);
fprintf(fout, '[END HEADER]\n\n');

% write edge map
fprintf(fout, '[BEGIN EDGEMAP]\n');
fprintf(fout, '# Format :  [EID] [Sub_Pixel_Pos] Sub_Pixel_Dir Strength\n');
for j=1:size(edges,1)
  % when saving edges, we make posions -1, because in cxx, index start
  % from 0, but in matlab index starts from 1
  s = ['[' num2str(j-1) ']' '[' num2str(edges(j,1)) ', ' num2str(edges(j,2)) ...
      ']  ' num2str(edges(j,3)),'   ',num2str(edges(j,4))];
  fprintf(fout, '%s\n', s);
end
fprintf(fout, '[END EDGEMAP]\n\n');

fprintf(fout, '[BEGIN CVLETMAP]\n');
fprintf(fout, ['# Format :  [EID] (number of curvelets)\n            ' ...
    '[edgels id] (forward) [ref_pt.x()] [ref_pt.y()] [ref_theta] [pt.x()] '...
    '[pt.y()] [theta k] length property\n']);
direction_str = 'BF';
for j = 1:size(edges,1)
    fprintf(fout, '<\n');
    % find out curvelets anchored at this edge
    cvlet_ind = find(chain(:,1)==j);
    cvlets = chain(cvlet_ind,2:end);
    fprintf(fout, ['[' num2str(j-1) ']']);
    fprintf(fout, ['(' num2str(size(cvlets,1)) ')\n']);
    for c = 1:size(cvlets,1)
        % if edge id = 0 it is invalid
        invalid_ind = find(cvlets(c,:)==0,1);
        if isempty(invalid_ind), invalid_ind = size(cvlets,2)+1; end
        format_str = repmat('%d ',1, invalid_ind-1);
        fprintf(fout, ['[' format_str '] (%s) [%f %f %f %f %f %f %f] %f %f\n'], ...
           cvlets(c,1:invalid_ind-1)-1, direction_str(cvletinfo(cvlet_ind(c),1)+1),...
           cvletinfo(cvlet_ind(c), 2:10));
    end
    fprintf(fout, '>\n');
end
fprintf(fout, '[END CVLETMAP]\n');

fclose(fout);