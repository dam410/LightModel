function [d] = distance_plane_camera(N,S,h)
	d = -h-transpose(N)*S;
end
