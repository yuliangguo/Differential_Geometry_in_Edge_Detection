function [subpix_x, subpix_y, subpix_dir_x, subpix_dir_y, subpix_grad_x, subpix_grad_y] = NMS_token(Gx, Gy, G, mask, margin)
%function [subpix_x, subpix_y, subpix_dir_x, subpix_dir_y, subpix_grad_x, subpix_grad_y] = NMS_token(Gx, Gy, G, mask, margin)
%
% NMS: non-maximal suppression code. Also detects the subpixel location of
% the maxima by fitting a parabola.
%
% Parameters:
%   Gx, Gy: Vector in which to do the nms (can be different from the
%           gradient vector if desired)
%   G:      The gradient magnitude surface (or any other surface)
%   margin: size of the border pixels to ignore
%
% Ouput:
%   subpix_x, subpix_y : subpixel location of the edgel token
%   subpix_dir_x, subpix_dir_y :  edgel orientation vector
%   subpix_grad_x, subpix_grad_y: gradient vector at edgel token
% 
%           6    7
%         -----------
%       5 |    |    | 8        
%         |    |    |
%         -----------
%       4 |    |    | 1
%         |    |    |
%         -----------
%           3     2
%
% (c) LEMS, Brown University
% Amir Tamrakar (amir_tamrakar@brown.edu)
% October 2007

%% Preprocess

%make sure Gmag, Gx and Gy are the same size

% output at least one edgel (bogus)
subpix_x(1)=0; 
subpix_y(1)=0;
subpix_dir_x(1)=0;
subpix_dir_y(1)=0;
subpix_grad_x(1)=0;
subpix_grad_y(1)=0;

%% Perform non-max suppression at every point
[size_y size_x] = size(G);
i = 1; %edgel token counter

for x=margin+2:size_x-(margin+2)
    for y=margin+2:size_y-(margin+2)
        gx = Gx(y,x);
        gy = Gy(y,x);

        if(mask(y,x)==0)
            continue;
        end

        if(abs(gx) < 10e-6 && abs(gy) < 10e-6) % invalid direction
            continue;
        end

        if(gx >= 0 && gy >= 0) %first quadrant
            if(gx >= gy)
                face = 1;
            else
                face = 2;
            end
        elseif(gx < 0 && gy >= 0)
            if(abs(gx) < gy)
                face = 3;
            else
                face = 4;
            end
        elseif(gx < 0 && gy < 0)
            if(abs(gx) >= abs(gy))
                face = 5;
            else
                face = 6;
            end
        elseif(gx >= 0 && gy < 0)
            if(gx < abs(gy))
                face = 7;
            else
                face = 8;
            end
        end

        % determine the three values for deciding whether this is a
        % max pixel and also to fit the parabola
        dir = [gx gy];
        dir = dir ./ norm(dir);
        f = G(y,x);
   
        if(face == 1)
            d1 = dir(2)/dir(1);
            fp = G(y,x+1) * (1-d1) + G(y+1,x+1) * d1;
            fm = G(y,x-1) * (1-d1) + G(y-1,x-1) * d1;
        elseif(face == 2)
            d1 = dir(1)/dir(2);
            fp = G(y+1,x) * (1-d1) + G(y+1,x+1) * d1;
            fm = G(y-1,x) * (1-d1) + G(y-1,x-1) * d1;
        elseif(face == 3)
            d1 = -dir(1)/dir(2);
            fp = G(y+1,x) * (1-d1) + G(y+1,x-1) * d1;
            fm = G(y-1,x) * (1-d1) + G(y-1,x+1) * d1;
        elseif(face == 4)
            d1 = -dir(2)/dir(1);
            fp = G(y,x-1) * (1-d1) + G(y+1,x-1) * d1;
            fm = G(y,x+1) * (1-d1) + G(y-1,x+1) * d1;
        elseif(face == 5)
            d1 = dir(2)/dir(1);
            fp = G(y,x-1) * (1-d1) + G(y-1,x-1) * d1;
            fm = G(y,x+1) * (1-d1) + G(y+1,x+1) * d1;
        elseif(face == 6)
            d1 = dir(1)/dir(2);
            fp = G(y-1,x) * (1-d1) + G(y-1,x-1) * d1;
            fm = G(y+1,x) * (1-d1) + G(y+1,x+1) * d1;
        elseif(face == 7)
            d1 = -dir(1)/dir(2);
            fp = G(y-1,x) * (1-d1) + G(y-1,x+1) * d1;
            fm = G(y+1,x) * (1-d1) + G(y+1,x-1) * d1;
        elseif(face == 8)
            d1 = -dir(2)/dir(1);
            fp = G(y,x+1) * (1-d1) + G(y-1,x+1) * d1;
            fm = G(y,x-1) * (1-d1) + G(y+1,x-1) * d1;
        end
        
        s = sqrt(1+d1^2);

        % max test 
        %   one can remove this test altogether and rely on the
        %   location of the max to decide whether to make an edgel token
        if((G(y,x) >  fm && G(y,x) >  fp) || ... %abs max
           (G(y,x) >  fm && G(y,x) >= fp) || ... %relaxed max
           (G(y,x) >= fm && G(y,x) >  fp))

            % fit parabola
            A = (fm+fp-2*f)/(2*s^2);
            B = (fp-fm)/(2*s);
            C = f;

            s_star = -B/(2*A); %location of max
            max_f = A*s_star*s_star + B*s_star + C; % value of max
  
            if(abs(s_star) <= sqrt(2)) % significant max is within a pixel
                % store the edgel token :
                % store the subpixel location
                subpix_x(i) = x + s_star * dir(1);
                subpix_y(i) = y + s_star * dir(2);

                %store the orientation of this edgel
                subpix_dir_x(i) = -dir(2);
                subpix_dir_y(i) = dir(1);

                %output gradient vector at max
                subpix_grad_x(i) = max_f*dir(1);
                subpix_grad_y(i) = max_f*dir(2);

                i = i+1;
            end %valid token
        end % max located

    end %y
end %x
