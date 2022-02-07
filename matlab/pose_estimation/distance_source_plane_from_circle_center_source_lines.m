function [h,d] = distance_source_plane_from_circle_center_source_lines(K_im,Xc,N,S)
	X_ray = inv(K_im)*[Xc;1];
	L_cell = lines_2p_to_plucker([zeros(3,1),S],[X_ray,N+S]);
	[Xs,err] = closest_point_lines(L_cell);
	h = dot(N,Xs) - dot(N,S);
	d = -dot(Xs,N);
end
