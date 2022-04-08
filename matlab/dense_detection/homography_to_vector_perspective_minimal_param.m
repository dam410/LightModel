% Initialize the perspectivity using affine considering a fronto parallel
%	view of the surface
function [h_vect_perspectivity,norm_t] = vector_init_to_perspectivity(h_vect_init,invK)
	% Calculate the inverse of the homography in normalized coordinates
	H = inv(vector_to_homography_perspective(h_vect_init)*inv(invK));
	N = cross(H(:,1),H(:,2));
	normN = norm(N);
	t = H(:,3)/norm(N);
	norm_t = norm(t);
	t_normalized = t./norm_t;
	N = N/norm(N);
	[alpha,beta] = normal_to_angle(N);
	h_vect_perspectivity = [alpha;beta;t_normalized(1);t_normalized(2);t_normalized(3)];
	vector_to_homography_perspective(h_vect_init);
	vector_to_perspectivity(h_vect_perspectivity,invK);
end
