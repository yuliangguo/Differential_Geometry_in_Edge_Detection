% Third-order Edge Detection Toolbox by Amir Tamrakar
%
% (c) LEMS, Brown University
% Amir Tamrakar (amir_tamrakar@brown.edu)
% October 2007
% 
% Various Resample Operators
%   G_2d_op -    Resampled Gaussian Filter  
%   Gx_2d_op -   Resampled Gx Filter
%   Gy_2d_op -   Resampled Gy Filter
%   Gxx_2d_op -  Resampled Gxx Filter
%   Gxy_2d_op -  Resampled Gxy Filter
%   Gyy_2d_op -  Resampled Gyy Filter
%   Gxxy_2d_op - Resampled Gxxy Filter
%   Gxyy_2d_op - Resampled Gxyy Filter
%   Gxxx_2d_op - Resampled Gxxx Filter
%   Gyyy_2d_op - Resampled Gyyy Filter
%   
% Resampled Convolution 
%   resampled_filter_2d - Interpolation filter using resampled filters
%
% Edge Detection
%   third_order_edge_detector - compute third-order edges
%   NMS_token - non-max suppression to produce edgel tokens
%
% Utility 
%   disp_edg - display the edgemap with line segments
%   example_usage - shows how to use the edge detector
