% Transform a plane pose into + ray direction of the circle center into an homography
%	h_vect = [alpha,beta,theta,gamma]
function [H] = vector_to_perspectivity_minimal_param(h_vect,invK)
	N = angle_to_normal(h_vect(1),h_vect(2));
	t = angle_to_normal(h_vect(3),h_vect(4));
	base_R = null(N*transpose(N));
	%if dot(cross(base_R(:,1),base_R(:,2)),N)<0
	R = [base_R(:,[1,2]),N];
	%else
	%	R = [base_R(:,[1,2]),N];
	%end
	M = [R(:,1:2),t];
	H = inv([R(:,1:2),t])*invK;
end
