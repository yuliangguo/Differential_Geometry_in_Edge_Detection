function save_edg(filename, edg, dim)
%save an .edg file

fid = fopen(filename, 'w');

% output header
fprintf(fid,'# EDGE_MAP v3.0\n');
fprintf(fid,'\n');
fprintf(fid,'# Format :  [Pixel_Pos]  Pixel_Dir Pixel_Conf  [Sub_Pixel_Pos] Sub_Pixel_Dir Sub_Pixel_Conf Sub_Pixel_Conf \n');
fprintf(fid,'\n');

%write out width and height info in the header
fprintf(fid, 'WIDTH=%d \n', dim(1));
fprintf(fid, 'HEIGHT=%d \n', dim(2));

edge_cnt = size(edg, 1);
fprintf(fid, 'EDGE_COUNT=%d \n', edge_cnt);
fprintf(fid,'\n\n');

for i=1:edge_cnt,
    % output each edge onto the file
    fprintf(fid, '[%d, %d]    %f %f  [%f, %f]  %f %f 0 \n', round(edg(i,1)), round(edg(i,2)), edg(i,3), edg(i,4), edg(i,1), edg(i,2), edg(i,3), edg(i,4));
end

%close the file
fclose(fid);