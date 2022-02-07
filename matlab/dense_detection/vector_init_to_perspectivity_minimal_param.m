% Initialize the perspectivity using affine considering a fronto parallel
%	view of the surface
function [h_vect_perspectivity,norm_t] = vector_init_to_perspectivity_minimal_param(h_vect_init,invK)
	% Calculate the decomposition of the homography
	%disp('Original init homography');
	%vector_to_homography_perspective(h_vect_init)
	[R,t] = decompose_homography(invK,vector_to_homography_perspective(h_vect_init));
	N = R(:,3);
	signN = sign(N(3));
	N = signN*N;
	t = t*signN;
	norm_t = norm(t);
	t_normalized = t./norm_t;
	[alpha,beta] = normal_to_angle(N);
	[gamma,delta] = normal_to_angle(t_normalized);
	h_vect_perspectivity = [alpha;beta;gamma;delta];
	%%inv_H_estimated = invK*inv(H_estimated);
	%%%inv_H_estimated = inv_H_estimated/norm(inv_H_estimated(:,1))
	%%%norm(inv_H_estimated(:,1))
	%%%norm(inv_H_estimated(:,2))
	%disp('Closest perspectivity');
	%H_estimated = vector_to_perspectivity_minimal_param(h_vect_perspectivity,invK)
end
