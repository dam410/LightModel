% Function that evaluate the photometric difference between en estimated photometric
%	model based on plane perspective projection of a concentric photometric model
%	following a spline model parametered with control points (I_pt,r_vect)
% Input parameters :
%	I : Complete image (not in used)
%	pt_in : Points coordinates (in image) that are being interpolated
%	I_pt : Points intensity value in pt_in
%	I_vect : Intensity values used for calculating the spline function
%	H : Homography representing the perpsective projection of the plane
%		Note that it can be affine
%	r_vect : Square radius used to control the spline
% Output parameters :
%	err : Difference between I_pt and model interpolated values
function [err] = evaluate_error_dense_detect(I,pt_in,I_pt,I_vect,H,r_vect)
        nb_pt = size(pt_in,1);
        % Project the points
        proj_pt = [pt_in,ones(nb_pt,1)]*transpose(H);
        proj_pt = proj_pt./proj_pt(:,3);
        r_pt = proj_pt(:,1).^2+proj_pt(:,2).^2;
        I_pt_proj = spline(r_vect,I_vect,r_pt);
        err = I_pt_proj-I_pt;
end

