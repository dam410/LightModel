function [pts_cell,I_pt_cell,I_vect_cell] = refinement_format(data,I,polys_2D,options)
	arguments
		data (1,1) struct
		I double
		polys_2D cell
		options.GlobalControlPoints (1,1) logical = false
		options.NB_R (1,1) double = 5
	end
	n_p = length(polys_2D);
	min_I_vect = Inf;
	max_I_vect = -Inf;
	I_pt_cell = cell(1,n_p);
	pts_cell = cell(1,n_p);
	I_vect_cell = cell(1,n_p);
	%figure, imshow(I/255);
	%hold on;
	for i_p = 1:n_p
		[pt_in_poly,I_pt,pt_barycenter,I_vect,ind_pt] = ...
			img_points_from_poly(I,polys_2D{i_p},options.NB_R);
		I_vect_cell{i_p} = I_vect;
		min_I_vect = min(min_I_vect,min(I_pt(:)));
		max_I_vect = max(max_I_vect,max(I_pt(:)));
		I_pt_cell{i_p} = I_pt;
		pts_cell{i_p} = pt_in_poly;
		%figure, plot3(pts_cell{i_p}(:,1),pts_cell{i_p}(:,2),I_pt_cell{i_p},'+r');
		%plot(polys_2D{i_p}(1,:),polys_2D{i_p}(2,:),'-r');
	end
	if options.GlobalControlPoints
		I_vect_cell = {min_I_vect:+(max_I_vect-min_I_vect)/(options.NB_R-1):max_I_vect};
	end
end
