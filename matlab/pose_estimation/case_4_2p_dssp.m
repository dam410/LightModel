function [dssp1,dssp2,psp1,psp2,ps] = case_4_2p_dssp(data)

	% Known parameters
	dssp1 = data.groundtruth.ScenePlaneDistanceSource{1};
	dssp2 = data.groundtruth.ScenePlaneDistanceSource{2};
	% Estimation
	[Xc1,N1] = plane_orientation_from_circular_contours(data.K,data.T_cam,...
		data.isocontour.CurveParameters{1});
	L_1 = source_line_from_distance_source_plane(data.K,Xc1,N1,dssp1);
	[Xc2,N2] = plane_orientation_from_circular_contours(data.K,data.T_cam,...
		data.isocontour.CurveParameters{2});
	L_2 = source_line_from_distance_source_plane(data.K,Xc2,N2,dssp2);
	% Source is calculated
	[ps,err] = closest_point_lines({L_1,L_2});
	d1 = distance_plane_camera(N1,ps,dssp1);
	d2 = distance_plane_camera(N2,ps,dssp2);
	psp1 = [N1;d1];
	psp2 = [N2;d2];
end
