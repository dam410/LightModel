% Calculate the normal direction from the circle center to the source
function [SX] = dir_S_X(K,S,Xc,N,d)
	X = circle_center_position(K,Xc,N,d);
	SX = (S-X)/norm(S-X);
end
