% Calculate virtual points from the homography which fit the homography
%	if image points are given, we use them instead
%function [points_iso,curves_iso] = isocontours_from_homography(poly_2D_mask,pt_in_poly,H_opt,N_circ)
function [points_iso,curves_iso] = isocontours_from_homography(poly_2D,H_opt,N_circ,I_pt,pt_im,spline_vector)
	% Just a trick to use only 3 arguments
	if nargin < 4
	end
	%[r_circle] = effective_radius(pt_in_poly,H_opt,N_circ);
	[r_circle] = effective_radius(transpose(poly_2D),H_opt,N_circ);
	t = 0:(2*pi/5000):2*pi;
	inv_H_opt = inv(H_opt);
	points_iso = cell(1,N_circ);
	curves_iso = cell(1,N_circ);
	for i_circ = 1:length(r_circle)
		% First calculate the conic 
		E = transpose(H_opt)*diag([1,1,-r_circle(i_circ)^2])*H_opt;
		% Calculate its parametric representation
		param_ell = ellipse2param(E);
		Pt_circ = inv_H_opt*[r_circle(i_circ)*cos(t);...
			r_circle(i_circ)*sin(t);ones(1,length(t))];
		Pt_circ = Pt_circ./Pt_circ(3,:);
		% Calculate the points
		%[ellipse_param,cov_ell_param] = ellipseFromPoints(transpose(Pt_circ([2,1],:)));
		curves_iso{i_circ} = param_ell;
		% Calculate the points which lie inside the polygons
		ind_Pt_circ = inpolygon(Pt_circ(1,:),Pt_circ(2,:),...
			poly_2D(1,:),poly_2D(2,:));
		points_iso{i_circ} = transpose(Pt_circ([2,1],ind_Pt_circ));
	end
end
