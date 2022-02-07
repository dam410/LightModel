function [h] = distance_source_plane_from_circle_center_source(K_im,Xc,N,S)
	% S = l.X - h.N
	if norm(S) < 1e-6
		% We choose a h so that the point is placed at 1 unit
		%	from the camera
		h = 0;
	else
		X_ray = inv(K_im)*[Xc;1];
		cross_X = cross_antisym(X_ray);
		XS = cross_X*S;
		XN = cross_X*N;
		h = -norm(XS)/norm(XN)*sign(dot(XS,XN));
	end
end
