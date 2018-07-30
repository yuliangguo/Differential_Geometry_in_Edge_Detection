function [edg, edgemap, thetamap] = load_edg(filename, opt)
%load an .edg file
% edg is Matlab coodinates

if nargin<2, opt=1; end

fid = fopen(filename);

ver = 1;
edge_cnt = 1;
width = 0;
height=0;

%scan the file to read the edge information
while 1
    lineBuffer = fgetl(fid);
    if ~ischar(lineBuffer), break, end

    if (strncmp(lineBuffer, '# EDGE_MAP v2.0', length('# EDGE_MAP v2.0')))
        ver = 2;
        continue;
    end
 
    if (strncmp(lineBuffer, '# EDGE_MAP v3.0', length('# EDGE_MAP v3.0')))
        ver = 2;
        continue;
    end
    %ignore other comment lines and empty lines
    if (length(lineBuffer)<2 | lineBuffer(1)=='#')
        continue;
    end

    if (strncmp(lineBuffer, 'WIDTH=', length('WIDTH=')))
        width = strread(lineBuffer,'WIDTH=%d');
        continue;
    elseif (strncmp(lineBuffer, ' WIDTH=', length(' WIDTH=')))
        width = strread(lineBuffer,' WIDTH=%d');
        continue;
    end
    
    if (strncmp(lineBuffer, 'HEIGHT=', length('HEIGHT=')))
        height = strread(lineBuffer,'HEIGHT=%d');
        continue;
    elseif (strncmp(lineBuffer, ' HEIGHT=', length(' HEIGHT=')))
        height = strread(lineBuffer,' HEIGHT=%d');
        continue;
    end

    % read the line with the edge count info
    if (strncmp(lineBuffer, 'EDGE_COUNT=', length('EDGE_COUNT=')))
        num_edges = strread(lineBuffer,'EDGE_COUNT=%d');
        if (ver==1)
            edg = zeros(num_edges, 4);
        elseif (ver==2)
            edg = zeros(num_edges, 5);
        end
        edgemap = zeros(height, width);
        thetamap = zeros(height, width);
        continue;
    elseif (strncmp(lineBuffer, ' EDGE_COUNT=', length(' EDGE_COUNT=')))
        num_edges = strread(lineBuffer,' EDGE_COUNT=%d');
        if (ver==1)
            edg = zeros(num_edges, 4);
        elseif (ver==2)
            edg = zeros(num_edges, 5);
        end
        edgemap = zeros(height, width);
        thetamap = zeros(height, width);
        continue;
    end

    % the rest should have data
    if (ver==1)
        % there are two variations of this file in existence
        if (strncmp(lineBuffer, 'EDGE : ', length('EDGE : ')))
            [ix, iy, idir, iconf, x, y, dir, conf] = strread(lineBuffer, 'EDGE :  [%d, %d]    %f %f   [%f, %f]   %f %f');
            if(conf==0)
                conf = 2;
            end
            edg(edge_cnt,:) = [x y dir conf];
            edgemap(round(y+1), round(x+1)) = conf;
            thetamap(round(y+1), round(x+1)) = dir;
            edge_cnt = edge_cnt + 1;
        else
            [ix, iy, idir, iconf, x, y, dir, conf] = strread(lineBuffer,' [%d, %d]   %f %f  [%f, %f]  %f %f');
            if(conf==0)
                conf = 2;
            end
            edg(edge_cnt,:) = [x y dir conf];
            edgemap(round(y+1), round(x+1)) = conf;
            thetamap(round(y+1), round(x+1)) = dir;
            edge_cnt = edge_cnt + 1;
        end
    elseif (ver==2)
        [ix, iy, idir, iconf, x, y, dir, conf, d2f] = strread(lineBuffer,' [%d, %d]   %f %f  [%f, %f]  %f %f %f');
        if(isempty(d2f))
            d2f = 0;
        end
        edg(edge_cnt,:) = [x y dir conf d2f];
            if(conf==0)
                conf = 2;
            end
        if (opt==1),
            edgemap(round(y+1), round(x+1)) = conf;
            thetamap(round(y+1), round(x+1)) = dir;
        elseif (opt==2)
            edgemap(round(y+1), round(x+1)) = -d2f;
        end
        edge_cnt = edge_cnt + 1;
    end
end

edgemap = edgemap(1:height,1:width);
%close the file
fclose(fid);

% update edge coordinates to matlab coordinates
% edg(:,1:2) = edg(:,1:2)+1;