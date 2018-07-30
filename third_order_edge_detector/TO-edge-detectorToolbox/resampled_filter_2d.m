function [filt_sig, dx] = resampled_filter_2d(sig, filt_fun_2d, sigma, n)
% function [filt_sig, dx] = resampled_filter_2d(sig, filt_fun_2d, sigma, n)
%
% Filter the signal at subpixel locations using resampled filters.
%   Resampled fitlers are simply filters whose analytic centers have been
%   shifted.
%
%   The benefits of this operation are:
%     1) Increased sampling of the response surface allows more small scale
%        edges to be detected. NMS typically relies on 3 neighboring pixels
%        to make its decision and also fit a parabola to find the subpixel
%        location of the edge. Edges at at only 1 or 2-pixel wide do not
%        pass this test and are thus rejected. Subpixel sampling of the
%        reponses surface allows these to be retained.
%
%     2) Increased sampling also results in smoother response surfaces and
%        hence subpixel localization improves.
%
% Parameters:
%    sig:          original image to be filtered
%    filt_fun_2d:  a 2d filter kernel that can be resampled (see G_2d_op)
%    n:            2^n number of subdivions within the pixel 
%
% Output:
%   filt_sig : resulting response signal 
%              Note that filt_sig will be 2^(2n) times larger than sig.
%
%   dx:        density of interpolation  = 1/(2^n); 
%
% See also G_2d_op, Gx_2d_op etc
%
% (c) LEMS, Brown University
% Amir Tamrakar (amir_tamrakar@brown.edu)
% October 2007

dx = 1/(2^n);  %sampling density

% define array to hold the filtered sig
[r, c] = size(sig);
filt_sig = zeros(r*2^n+1, c*2^n+1);

% filter with shifted versions of the original filter
% we need 2^(2n) convolutions to completely compute it
for i=0:2^n-1,
    for j=0:2^n-1,
        %compute subpixel shifts
        dp_x = i/(2^n);
        dp_y = j/(2^n);
        
        % compute the filter at this shift 
        filt_resamp = filt_fun_2d(sigma, -dp_x, -dp_y);
        
        % filter sig with this filter
        filt_sig_p = imfilter(sig, filt_resamp, 'conv','replicate');
        
        % filter sig with 1-d versions
        %[filt_resamp, filt_x, filt_y] = filt_fun_2d(sigma, -dp_x, -dp_y);
        %filt_sig_px = imfilter(sig, filt_x, 'conv','replicate');
        %filt_sig_p2 = imfilter(filt_sig_px, filt_y', 'conv','replicate');
        
        filt_sig(j+1:2^n:r*2^n,i+1:2^n:c*2^n) = filt_sig_p;
    end
end

%chop off the excess at the boundaries
filt_sig = filt_sig(1:(r-1)*2^n+1, 1:(c-1)*2^n+1);