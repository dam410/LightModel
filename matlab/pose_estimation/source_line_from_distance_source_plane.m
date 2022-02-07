function [L_s] = source_line_from_distance_source_plane(K,Xc,N,h)
	X_line = -h*N;	
	d = inv(K)*[Xc;1];
	d = d/norm(d);
	m = cross(X_line,X_line+d);
	L_s = plucker_dm_to_matrix(d,m);
end
