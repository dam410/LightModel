% Calculate the geometric error between homographies
%	with invariance in rotation and scale
function [err,r_pt_1_normalized,r_pt_2_normalized] = homography_error(I,poly_2D,H_1,H_2)
	% Get only the point inside the polygon
	[pt_in] = img_points_from_poly(I,poly_2D,5);
	% Calculate their radius with both homographies
	[r_pt_1] = radius_from_H(pt_in,H_1);	
	[r_pt_2] = radius_from_H(pt_in,H_2);	
	% Normalized the radius
	r_pt_1_normalized = sqrt(r_pt_1)./mean(sqrt(r_pt_1));
	r_pt_2_normalized = sqrt(r_pt_2)./mean(sqrt(r_pt_2));
	err = r_pt_1_normalized-r_pt_2_normalized;
end
