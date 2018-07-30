function f = mTOdetector_features(img, sigma, radius, n, islogist)
% now only one channel: lumination

img_g = rgb2gray(img);
img_g = double(img_g);

% compute gaussian derivation response map for each scale
for s = 1:length(sigma)
    Ix(:,:,s)   = resampled_filter_2d(img_g, @Gx_2d_op, sigma(s), n);
    Iy(:,:,s)   = resampled_filter_2d(img_g, @Gy_2d_op, sigma(s), n);
    grad_mag(:,:,s) = sqrt(Ix(:,:,s).^2+Iy(:,:,s).^2);
end

[d1, d2, d3] = size(grad_mag);

% compute texture gradient maps
img = imresize(img, [d1 d2], 'cubic');
for i = 1:length(radius)
    tg_map = tg_features(img, radius(i));
    grad_mag(:,:,d3+i) = tg_map;
end

[~, ~, d3] = size(grad_mag);
grad_mag = reshape(grad_mag, [d1*d2, d3]);

f = grad_mag';
if(islogist)
    f= [ones(1, d1*d2); f];
end

