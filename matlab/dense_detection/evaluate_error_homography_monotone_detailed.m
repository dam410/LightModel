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
%	a_vect : Square radius used to control the spline
% Output parameters :
%	err : Difference between I_pt and model interpolated values
function [I_pt_proj,r_proj,err] = evaluate_error_homography_monotone_detailed(I,pt_in,I_pt,I_vect,H,a_vect)

	n_vect = length(a_vect);
	J_r_vect = tril(ones(n_vect));
        r_vect = J_r_vect*reshape(a_vect,n_vect,1);
        I_vect = reshape(I_vect,n_vect,1);

	% Calculate the squared radius
	[r_pt,J_r_pt] = radius_from_H(pt_in,H);
        nb_pt = length(r_pt);
	%% Calculate the spline
	%[pp,J_pp_mat] = monotone_spline(r_vect,I_vect);
	%% Calculate the intensity
	%I_pt_proj = ppval(pp,r_pt);
	r_proj = r_pt;
	try
		pp = monotone_spline(r_vect,I_vect);
		I_pt_proj = ppval(pp,r_pt);
	catch exception
		%disp('Problem of non increasing r_vect');
		I_pt_proj = spline(r_vect,I_vect,r_pt);
	end	
        err = I_pt_proj-I_pt;
end
