% Initialize the rotation using affine initialization (rough)
function [h_vect_rotation,scale] = vector_init_to_rotation_param(h_vect_init,invK)
	% Calculate the decomposition of the homography
	H_rough = invK*inv(vector_to_homography_perspective(h_vect_init));
	scale_mat = sqrt(H_rough(1,1)^2+H_rough(1,1)^2+H_rough(2,2)^2);
	H_rough_normalized = H_rough/scale_mat;
	N = H_rough_normalized(:,3);
	scale = norm(N);
	signN = sign(N(3));
	N = signN*N/norm(N);
	[alpha,beta] = normal_to_angle(N);
	h_vect_rotation = [alpha;beta];
end
