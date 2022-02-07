function [data] = add_isocontours_to_homogs_ambig(data,polys_2D,homogs_1,homogs_2,N_CIRC)
	nb_planes = length(polys_2D);
	points = cell(1,nb_planes);
	curves = cell(1,nb_planes);
	homographies = cell(1,nb_planes);
	for i_p = 1:nb_planes
		[points_iso_1,curves_iso_1] = isocontours_from_homography(...
			polys_2D{i_p},homogs_1{i_p},N_CIRC);
		[points_iso_2,curves_iso_2] = isocontours_from_homography(...
			polys_2D{i_p},homogs_2{i_p},N_CIRC);
		points{i_p} = {points_iso_1,points_iso_2};
		curves{i_p} = {curves_iso_1,curves_iso_2};
		homographies{i_p} = {homogs_1{i_p},homogs_2{i_p}};
	end
	data = data;
	isocontour = struct();
	isocontour.CurveParameters = curves;
	isocontour.Points = points;
	data.isocontour = isocontour;
	data.homography = homographies;
end
