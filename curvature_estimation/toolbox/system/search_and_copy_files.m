clc;
clear;

list_path = '~/Documents/MATLAB/Third_order/Data/BSDS300_data/images/train/';
src_path = '~/Documents/MATLAB/Third_order/Data/BSDS500_RES/train_curve/MEIM_5/';
dst_path = '~/Documents/MATLAB/Third_order/Data/BSDS300_RES/train_curve/MEIM_5/';
mkdir(dst_path);
ext = '.cem';

files = dir([list_path]);
files = files(3:end);

count = 0;
for i=1:length(files)
   fprintf(1,'file: %d/%d \n', i, length(files));
   targetFile = [src_path files(i).name(1:end-4) ext];
   if exist(targetFile, 'file')
       count = count+1;
       copyfile(targetFile,[dst_path files(i).name(1:end-4) ext]);
   else
       fprintf(2,'\t file: %s not found\n', [files(i).name(1:end-4) ext]);
   end
end

fprintf(1, '%d/%d files found and copied\n', count, length(files));