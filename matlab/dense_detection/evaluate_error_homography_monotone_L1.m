% Function that evaluate the photometric difference between en estimated photometric
%	model based on plane perspective projection of a concentric photometric model
%	following a spline model parametered with control points (I_pt,r_vect)
%	with a L1 norm !!
% Input parameters :
%	I : Complete image (not in used)
%	pt_in : Points coordinates (in image) that are being interpolated
%	I_pt : Points intensity value in pt_in
%	I_vect : Intensity values used for calculating the spline function
%	H : Homography representing the perpsective projection of the plane
%		Note that it can be affine
%	a_vect : Square radius used to control the spline
% Output parameters :
%	norm_L1 : Norm L1 of the difference between I_pt and model interpolated values
%	gradient_norm_L1 : Gradient of the L1 norm of the error with variable parameters, [H(:);r_vect(:)]
function [norm_L1,gradient_norm_L1] = evaluate_error_homography_monotone_L1(I,pt_in,I_pt,I_vect,H,a_vect)
	n_pt = length(I_pt);
	n_vect = length(a_vect);
	if nargout > 1
		[err,J_err_full] = evaluate_error_homography_monotone(I,pt_in,I_pt,I_vect,H,a_vect);
		[norm_L1,gradient_norm_L1] = error_norm_L1(err,J_err_full);
	else
		[err] = evaluate_error_homography_monotone(I,pt_in,I_pt,I_vect,H,a_vect);
		[norm_L1] = error_norm_L1(err);
	end 
end

