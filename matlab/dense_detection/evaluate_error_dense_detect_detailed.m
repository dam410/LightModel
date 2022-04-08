% See function evaluate_error_dense_detect, does the same but add more output detailed
% Output parameters:
%	I_pt_proj : Intensity values interpolated with the model
%	r_proj : control square radius for the spline
%	err : Same as function
function [I_pt_proj,r_proj,err] = evaluate_error_dense_detect_detailed(I,pt_in,I_pt,I_vect,H,r_vect)
        nb_pt = size(pt_in,1);
        % Project the points
        proj_pt = [pt_in,ones(nb_pt,1)]*transpose(H);
        proj_pt = proj_pt./proj_pt(:,3);
        r_proj = proj_pt(:,1).^2+proj_pt(:,2).^2;
        I_pt_proj = spline(r_vect,I_vect,r_proj);
        err = I_pt_proj-I_pt;
end

