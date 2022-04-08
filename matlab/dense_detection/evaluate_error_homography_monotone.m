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
function [err,J_err] = evaluate_error_homography_monotone(I,pt_in,I_pt,I_vect,H,a_vect)
	n_vect = length(a_vect);
	J_r_vect = tril(ones(n_vect));
	r_vect = J_r_vect*reshape(a_vect,n_vect,1);
	I_vect = reshape(I_vect,n_vect,1);
	% Calculate the squared radius
	[r_pt,J_r_pt] = radius_from_H(pt_in,H);
        [nb_pt,nb_param_h] = size(J_r_pt);
	%try
	if nargout > 1
		% With Jacobian calculation
		[pp,J_pp_mat,J_breaks] = monotone_spline(r_vect,I_vect);
		J_pp_mat_r = multiplication_tensor_matrix(J_pp_mat(:,:,1:n_vect),J_r_vect);
		J_breaks_r = J_breaks(:,1:n_vect)*J_r_vect;
		%% Calculate the intensity
		[I_pt_proj,J_err] = ppval_derivative_independent(pp,J_pp_mat_r,J_breaks_r,r_pt,J_r_pt);
		% Jacobian is built with parameters of r_vect first (but most of the
		%	time we want parameter of r_pt first)
		J_err = J_err(:,[n_vect+(1:nb_param_h),1:n_vect]);
	else
		% Without jacobian calculation
		pp = monotone_spline(r_vect,I_vect);
		I_pt_proj = ppval(pp,r_pt);
	end
	%%catch exception
	%	%disp('Problem of non increasing r_vect');
	%	I_pt_proj = spline(r_vect,I_vect,r_pt);
	%	J_err = zeros(length(I_pt_proj),nb_param_h+n_vect);
	%end
	% Calculate the final error and its jacobian
        err = I_pt_proj-I_pt;
end

