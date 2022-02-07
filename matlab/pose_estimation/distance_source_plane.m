function [h] = distance_source_plane(N,S,d)
	h = -(d+transpose(N)*S);
end
