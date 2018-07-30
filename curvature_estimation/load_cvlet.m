function [ chain, config, edg, opts] = load_cvlet(filename)
%load a .cem file (v2.0)
% In cem, coodinates keep the original
% ref_pt.x() ref_pt.y() are enforced to save K_min and K_max !!!!

fid = fopen(filename);

if (fid<0)
    fprintf(2,'File not found.');
    chain = [];
    config = [];
    edg = [];
    opts = strct();
    return;
end


lineBuffer = fgetl(fid); %read in the first line

if (~strncmp(lineBuffer, 'CVLET_MAP v2.0', length('CVLET_MAP v2.0')))
    return; %wrong version
end

% scan the file to read the contour information
% read the header info
while 1,
    lineBuffer = fgetl(fid);
    if (strncmp(lineBuffer, '[END HEADER]', length('[END HEADER]'))), break, end

    %ignore comment lines and empty lines
    if (length(lineBuffer)<2 | lineBuffer(1)=='#')
        continue;
    end

    % read the header info
    % width
    if (strncmp(lineBuffer, 'WIDTH=', length('WIDTH=')))
        opts.w = strread(lineBuffer,'WIDTH=%d'); %#ok<*DSTRRD>
        continue;
    end
    % height
    if (strncmp(lineBuffer, 'HEIGHT=', length('HEIGHT=')))
        opts.h = strread(lineBuffer,'HEIGHT=%d');
        continue;
    end
    % # of edges
    if (strncmp(lineBuffer, 'EDGE_COUNT=', length('EDGE_COUNT=')))
        edgecount = strread(lineBuffer,'EDGE_COUNT=%d');
        continue;
    end 
    % curve model type
    if (strncmp(lineBuffer, 'CM TYPE=', length('CM TYPE=')))
        opts.cmtype = strread(lineBuffer,'CM TYPE=%s');
        opts.cmtype = opts.cmtype{1};
        continue;
    end
    % radius of neighborhood
    if (strncmp(lineBuffer, 'N_RADIUS=', length('N_RADIUS=')))
        opts.nrad = strread(lineBuffer,'N_RADIUS=%f');
        continue;
    end
    % position uncertainty
    if (strncmp(lineBuffer, 'DX=', length('DX=')))
        opts.dx = strread(lineBuffer,'DX=%f');
        continue;
    end
    % orientation uncertainty
    if (strncmp(lineBuffer, 'DT=', length('DT=')))
        opts.dt = strread(lineBuffer,'DT=%f');
        continue;
    end
    % token length
    if (strncmp(lineBuffer, 'TOKEN_LENGTH=', length('TOKEN_LENGTH=')))
        opts.token_len = strread(lineBuffer,'TOKEN_LENGTH=%f');
        continue;
    end
    % max size to group
    if (strncmp(lineBuffer, 'MAX_SIZE_TO_GROUP=', length('MAX_SIZE_TO_GROUP=')))
        opts.max_size_to_group = strread(lineBuffer,'MAX_SIZE_TO_GROUP=%d');
        continue;
    end
end

if ~isfield(opts,'max_size_to_group')
    warning('Maxium size to group does not recorded. Set as 15.')
    opts.max_size_to_group = 15;
end

% initialize the data structure
edg = zeros(edgecount,4);
chain = cell(edgecount,1);
config = cell(edgecount,1);
max_cvlet_len = 0;

while 1,
    lineBuffer = fgetl(fid);
    if ~ischar(lineBuffer)
        fclose(fid);
        return;
    end

    %ignore comment lines and empty lines
    if (length(lineBuffer)<2 | lineBuffer(1)=='#')
        continue;
    end
    %read the edgemap block
    if (strncmp(lineBuffer, '[BEGIN EDGEMAP]', length('[BEGIN EDGEMAP]')))
        %read the next line with comments
        lineBuffer = fgetl(fid);
        %read in all the edges
        for i=1:edgecount
            lineBuffer = fgetl(fid);
            [~, x, y, dir, conf] = strread(lineBuffer,'[%d] [%f, %f]   %f %f');
            edg(i,:) = [x, y, dir, conf];
        end
        continue;
    end

    % read the curvelet block
    if (strncmp(lineBuffer, '[BEGIN CVLETMAP]', length('[BEGIN CVLETMAP]')))
        while(1)
            % read in all the curvelets
            lineBuffer = fgetl(fid);
            if ~ischar(lineBuffer)
                config = cat(1,config{:});
                chain = cat(1,chain{:});
                if opts.max_size_to_group > max_cvlet_len
                    chain = chain(:,1:max_cvlet_len+1);
                    opts.max_size_to_group = max_cvlet_len;
                end
                fclose(fid);
                return
            end
            % read curvelets of an edgel
            if(strncmp(lineBuffer, '<', length('<')))
                lineBuffer = fgetl(fid);
                [id,num] = strread(lineBuffer,'[%d] (%d)');
                % c style id starts from 0
                id = id+1;
                % for edge ID
                chain{id} = zeros(num,opts.max_size_to_group+1);
                chain{id}(:,1) = id; 
                % for configurations
                config{id} = zeros(num,10);
                for i = 1:num
                    %read in the list of edge ids
                    lineBuffer = fgetl(fid);
                    strcell = strsplit(lineBuffer,{'[',']'});
                    edgeID = strread(strcell{2},'%d')';
                    % c style id starts from 0
                    edgeID = edgeID+1;
                    % ensure the direction of curvelet is forward. If not, reverse
                    % the edge ID chain.
%                     if(strcell{3}(3)~='F')
%                         edgeID = fliplr(edgeID);
%                     end
                    chain{id}(i,2:numel(edgeID)+1) = edgeID;
                    max_cvlet_len = max([max_cvlet_len, numel(edgeID)]);
                    
                    % configuration of a curvelet: [isForward ref_pt.x()] [ref_pt.y()] 
                    % [ref_theta] [pt.x()] [pt.y()] [theta] k]
                    if strcell{3}(3)=='F', config{id}(i,1) = 1; end
                    config{id}(i,2:8) = strread(strcell{4},'%f','delimiter',' ');
                    config{id}(i,9:10) = strread(strcell{5},'%f','delimiter',' ');
                end
            end
        end
    end
end





