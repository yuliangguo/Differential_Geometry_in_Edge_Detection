function clearFolder(src_path)
% clear the temp file start with '._' in current folder
% and its subfolders
%
% USAGE
%  clearFolder(src_path, isSubfolder)
%
% INPUTS
%  src_path   - [string] the path of the folder to be cleared
%
% Code written by Xiaoyan Li, 2015.

files = dir(src_path);
files = files(3:end);

if any([files.isdir])
    subFolder = {files([files.isdir]).name};
    for j = 1:length(subFolder)
        subFolder{j} = [subFolder{j} '/'];
    end
    subFolder = {subFolder{:}, ''};
else
    subFolder = {''};
end

k = 0;

for j = 1:length(subFolder)
    tempinfo = dir([src_path subFolder{j} '._*']);
    for i = 1:length(tempinfo)
        delete([src_path subFolder{j} tempinfo(i).name])
        k = k+1;
    end
end
disp(['remove ' num2str(k) ' temp files']);

end
