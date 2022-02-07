% This function use the photometrically optimal homographies to calculate 
%	discrete isocontours in the image
%	that can be later use to estimate the pose of the planes and PLS
%	or as points correspondences for global geometric refinement
function [data] = add_isocontours_to_homog(data,polys_2D,homogs,N_CIRC)
	nb_planes = length(polys_2D);
	points = cell(1,nb_planes);
	curves = cell(1,nb_planes);
	homographies = cell(1,nb_planes);
	for i_p = 1:nb_planes
		[points_iso,curves_iso] = isocontours_from_homography(...
			polys_2D{i_p},homogs{i_p},N_CIRC);
		points{i_p} = {points_iso};
		curves{i_p} = {curves_iso};
		homographies{i_p} = {homogs{i_p}};
	end
	data = data;
	isocontour = struct();
	isocontour.CurveParameters = curves;
	isocontour.Points = points;
	data.isocontour = isocontour;
	data.homography = homographies;
end
