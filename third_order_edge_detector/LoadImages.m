function [imgpath] = LoadImages(path,ext)
str = path;
pt = dir(str);
foldname = []; 
foldname{1}='';
k = 1;
imgpath = {};
for i = 1 : length(pt)
    if strcmp(pt(i).name, '.') | strcmp(pt(i).name, '..')
        % compares the strings and returns logical 1(true) if they are identical  
        continue;
    else
        k = k + 1;
        foldname{k} = pt(i).name;
    end
end
for i = 1 : length(foldname)
    temp = strcat(str, foldname{i}, '/*', ext);
    temp1 = dir(temp);
    for j = 1 : length(temp1)
        imgpath{i, j} = strcat(str, foldname{i}, '/', temp1(j).name);
    end
end


