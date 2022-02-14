function [L_s] = source_line_from_plane_pose(K,Xc,N,d)
	% Calculate the position of the circle center coordinates
	X_ray = inv(K)*[Xc;1];
	lambda = -d/(transpose(N)*X_ray);
	X_ray = lambda*X_ray;
	d = N;
	m = cross(X_ray,X_ray+d);
	L_s = plucker_dm_to_matrix(d,m);
	%L_s = [0,-d(1),-d(2),-d(3);d(1),0,-m(1),-m(2);d(2),m(1),0,-m(3);d(3),m(2),m(3),0];
end
