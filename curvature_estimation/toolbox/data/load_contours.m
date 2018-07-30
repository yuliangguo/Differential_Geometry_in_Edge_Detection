function cem = load_contours(filename)
%load a .cem file (v2.0)

fid = fopen(filename);

if (fid<0)
    disp('File not found.');
    cem = {};
    return;
end

%cem structure
cem = cell(3,1);

% cells to store the contour information
contours = {};
num_contours = 0;

lineBuffer = fgetl(fid); %read in the first line

if (~strncmp(lineBuffer, '.CEM v2.0', length('.CEM v2.0')))
    return; %wrong version
end

%scan the file to read the contour information
while 1,
    lineBuffer = fgetl(fid);
    if ~ischar(lineBuffer), break, end

    %ignore comment lines and empty lines
    if (length(lineBuffer)<2 | lineBuffer(1)=='#')
        continue;
    end

    % read the line with the edgemap size info
    if (strncmp(lineBuffer, 'size=', length('size=')))
        [w, h] = strread(lineBuffer,'size=[%d %d]');
        cem{1} = [w h];
        continue;
    end

    %read the edgemap block
    if (strncmp(lineBuffer, '[Edgemap]', length('[Edgemap]')))
        %read the next line with the edge count
        lineBuffer = fgetl(fid);
        num_edges = strread(lineBuffer,'count=%d');

        %read in all the edges
        for i=1:num_edges,
            lineBuffer = fgetl(fid);
            [x, y, dir, conf, d2f] = strread(lineBuffer,'(%f, %f)\t%f\t%f\t%f');
            edg(i,:) = [x, y, dir, conf, d2f];
        end
        continue;
    end

    % read the contours block
    if (strncmp(lineBuffer, '[Contours]', length('[Contours]')))
        %read the next line with the contour count
        lineBuffer = fgetl(fid);
        num_contours = strread(lineBuffer,'count=%d');

        %read in all the contours
        for i=1:num_contours,
            %read in the list of edge ids
            lineBuffer = fgetl(fid);
            e_ids = strread(lineBuffer(2:length(lineBuffer)-1),'%d','delimiter',' ')';

            chain = [];
            for j=1:length(e_ids)
                chain(j,:) = [edg(e_ids(j)+1,:) e_ids(j)+1];
            end

            %store chain
            contours{i} = chain;
        end
        
        %save the contours
        cem{2} = contours;
        continue;
    end
    
    % read the contour properties block
    if (strncmp(lineBuffer, '[Contour Properties]', length('[Contour Properties]')))
        lineBuffer = fgetl(fid);
        
        %read in all the contour properties
        con_props = [];
        for i=1:num_contours,
            %read in the list of contour properties
            lineBuffer = fgetl(fid);
            
            %<len> <avg. str> <mean con> <Lstd> <Rstd> <avg. d2f> <avg. k> <max k>
            con_props(i,:) = strread(lineBuffer,'%f','delimiter',' ');
        end
        
        %save the properties
        cem{3} = con_props;
    end
end

%close the file
fclose(fid);

% %for debug
% colourmp = hsv(length(contours));    % HSV colour map with con_cnt entries
% colourmp = colourmp(randperm(length(contours)),:);  % Random permutation
% 
% figure;
% axis ij;
% for i = 1:length(contours)-1
%     line(contours{i}(:,1), contours{i}(:,2),'color',colourmp(i,:));
% end





