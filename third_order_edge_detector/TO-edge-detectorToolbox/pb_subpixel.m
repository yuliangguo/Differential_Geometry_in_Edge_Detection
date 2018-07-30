function edg = pb_subpixel(im, thresh, opt)
%
% compute subpixel edges from pb response
%
% opt: 0=pBG, 1:pTG, 2:pBGTG, 3:pCG, 4:pCGTG
%

% compute pb according to the option
switch(opt)
    case 0
        [pb, theta] = pbBG_soft(im);
    case 1
        [pb, theta] = pbTG_soft(im);
    case 2
        [pb, theta] = pbBGTG_soft(im);
    case 3
        [pb, theta] = pbCG_soft(im);
    case 4
%         [pb, theta] = pbCGTG(im); 
        [pb, theta] = pbCGTG_soft(im);
end

% Use NMS to detect the edgels - this NMS routine ouputs subpixel edgel tokens
margin = 4;
dirx = cos(theta+pi/2);
diry = sin(theta+pi/2);
[subpix_x, subpix_y, subpix_dir_x, subpix_dir_y, subpix_grad_x, subpix_grad_y] = NMS_token(dirx, diry,pb, pb>thresh, margin);
mag_e = sqrt(subpix_grad_x.^2+subpix_grad_y.^2);

% construct the edge map
edg = [subpix_x' subpix_y' atan2(subpix_dir_y',  subpix_dir_x')  mag_e'];