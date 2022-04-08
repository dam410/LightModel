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
%	J_err : Jacobian of the error with variable parameters, [H(:);r_vect(:)]
function [err,J_err] = evaluate_error_homography_monotone_proxy_only_x_spline(I,pt_in,I_pt,I_vect,H,as_vect)
	% Now we suppose that a_vect is only the sqrt of the full incremental function
	%	-> We need to update the jacobian accordingly
	[err,J_err_full] = evaluate_error_homography_monotone(I,pt_in,I_pt,I_vect,H,as_vect.^2);
	n_end = size(J_err_full,2);
	n_end_param = n_end-2*length(I_vect);
	n_end_xspline = n_end-length(I_vect);
	J_err_param = J_err_full(:,1:n_end_param);
	J_err_spline = J_err_full(:,(n_end_xspline+1):n_end);%.*transpose(as_vect);
	J_err = [J_err_param,J_err_spline];
end

