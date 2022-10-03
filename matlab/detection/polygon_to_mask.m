% Function to convert a polygon based on points into a contour in the image
%	* Suppose the polygon is convex !!
% 	* This takes into account the border of the image
%	* A parameter can be given to contract the shape to avoid
%	offset error which will make the active contour to stop on the edge of the shape
function [resulting_mask] = polygon_to_mask(poly_2D,I,contraction_offset)
	poly_2D
	EPS = 1e-10;
	size_I = size(I);
	[n,m] = size(I);
	nb_pt = size(poly_2D,2);
	%  Iterate over the points to define a new polygon which is inside the boundary
	%  of an image
	i_first_inbound = 0;
	inboundcheck = 0;
	while i_first_inbound<nb_pt & ~inboundcheck
		i_first_inbound = i_first_inbound+1;
		inboundcheck = in_bc(poly_2D(:,i_first_inbound),size_I);
	end
	if ~inboundcheck
		disp('No points to start, contour will be the border of the image');
		pts_contour = [1,size_I(1),size_I(1),1;...
			1,1,size_I(2),size_I(2)];
	end
	% Try to recreate a vector of indices starting from a point inside the image boundary
	if i_first_inbound==1
		index = [(i_first_inbound+1):(nb_pt)];
	else
		index = [(i_first_inbound+1):(nb_pt),1:(i_first_inbound-1)];
	end
	pt_from = poly_2D(:,i_first_inbound)
	pt_boundary = [pt_from];
	index
	for i = index
		[inboundcheck,pt_inter] = in_bc(poly_2D(:,i),size_I,pt_from);
		if inboundcheck
			pt_boundary = [pt_boundary,poly_2D(:,i)];
			pt_from = poly_2D(:,i);
		else
			pt_boundary = [pt_boundary,pt_inter];
		end
	end
	% Addind the first point once again at the end and the last at the beginning
	pt_boundary_cycle = [pt_boundary(:,end),pt_boundary,poly_2D(:,i_first_inbound)];
	pt_boundary_offset = [];
	% Smoothing by moving each points into the direction of its neighboughs
	for j = 2:(length(pt_boundary_cycle)-1)
		direction = (pt_boundary_cycle(:,j-1)-pt_boundary_cycle(:,j))+(pt_boundary_cycle(:,j+1)-pt_boundary_cycle(:,j));
		direction = direction/norm(direction);
		if norm(direction)<EPS
			direction = [0;0];
		end
		pt_boundary_offset = [pt_boundary_offset,contraction_offset*direction+pt_boundary_cycle(:,j)];
	end
	pt_boundary_offset = [pt_boundary_offset,pt_boundary_offset(:,1)];
	resulting_mask = poly2mask(pt_boundary_offset(1,:),pt_boundary_offset(2,:),m,n);
end


