% This function uses the known homography to calculate the radius with r_min and r_max set
%       to min and max values
% Adding the parameter I_pt to be able to setup r_vect a bit better
function [h_vect,r_vect] = initial_parameters_given_homography(I,pt_in_poly,I_vect,H)
        nb_pt = length(pt_in_poly);
        % Calculate the inverse to get the initial parameters for the homography
        invH = inv(H);
        invH = invH./invH(3,3);
        h_vect = transpose(invH(1:8));
        n_val = length(I_vect);
        % Project the points
        H = reshape([h_vect;1],3,3);
        proj_pt = [pt_in_poly,ones(nb_pt,1)]*transpose(H);
        proj_pt = proj_pt./proj_pt(:,3);
        r_proj = proj_pt(:,1).^2+proj_pt(:,2).^2;
        r_max = max(r_proj);
        r_min = min(r_proj);
        r_vect = transpose(r_min:(r_max-r_min)/(n_val-1):r_max);
end

