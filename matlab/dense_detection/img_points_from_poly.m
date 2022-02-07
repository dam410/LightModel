% From one image and the local polygon surface delimitation, we extract :
%	pt_in_poly : Points coordinate inside the polygon
%	I_pt : Their intensities values
%	pt_barycenter : Their barycenter
%	I_vect : Vectors of intensities ranging from min(I_pt) to max(I_pt) to interpolate
function [pt_in_poly,I_pt,pt_barycenter,I_vect,ind_pt] = img_points_from_poly(I,poly_2D,NB_R)
	
	[n,m] = size(I);

	% Use a mask to search only inside the polygon
	pt_barycenter = mean(poly_2D,2);
	direction_shrink = pt_barycenter-poly_2D;
	norm_direction_shrink = sqrt(sum(direction_shrink.^2));
	direction_shrink = direction_shrink./norm_direction_shrink;
	poly_2D_mask = poly_2D+10*direction_shrink;
	% Calculate the pixel points which lie inside the surface
	%	and the image (in min and max coordinates)
	min_x = max(1,min(m,floor(min(poly_2D_mask(1,:)))));
	max_x = max(1,min(m,ceil(max(poly_2D_mask(1,:)))));
	min_y = max(1,min(n,floor(min(poly_2D_mask(2,:)))));
	max_y = max(1,min(n,ceil(max(poly_2D_mask(2,:)))));
	[X,Y] = meshgrid(min_x:max_x,min_y:max_y);
	in_poly = inpolygon(X(:),Y(:),...
		poly_2D_mask(1,:),poly_2D_mask(2,:));
	pt_in_poly = [X(in_poly),Y(in_poly)];

	% Calculate the intensity values for selected pixels with min and max value 
	ind_pt = sub2ind([n,m],pt_in_poly(:,2),pt_in_poly(:,1));
	I_pt = I(ind_pt);
	[I_max,ind_max] = max(I_pt);
	[I_min,ind_min] = min(I_pt);

	% Select the intensity range
	I_vect = I_max:-(I_max-I_min)/(NB_R-1):I_min;

end
