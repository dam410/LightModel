function [P_s] = source_plane_from_circle_center_orientation(K,Xc,N)
	X_ray = inv(K)*[Xc;1];
	X_ray = X_ray/norm(X_ray);
	N_s = cross(X_ray,N);
	N_s = N_s/norm(N_s);
	P_s = [N_s;0];
end
