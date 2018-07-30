function s = LeastSquareObj(w, f, y)

% load('features_mTO_CRIM.mat');
% 
% f=f'; y=y';
% % normalize features to unit variance
% fstd = std(f);
% fstd = fstd + (fstd==0);
% f = f ./ repmat(fstd,size(f,1),1);

s = (f*w-y)'*(f*w-y);