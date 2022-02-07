function [S] = source_from_plane_pose_distance_source_plane(K,Xc,N,d,h)
	% Calculate the position of the circle center
	X = circle_center_position(K,Xc,N,d);
	S = X - h*N;
end
