function [dssp1,dssp2,psp1,psp2,ps] = case_6_1p_psp(data)

	% Known parameters
	psp1 = [data.groundtruth.ScenePlaneOrientation{1}(:,3);data.groundtruth.ScenePlanePosition{1}];
	psp2 = [data.groundtruth.ScenePlaneOrientation{2}(:,3);data.groundtruth.ScenePlanePosition{2}];
	% Estimation
	[Xc1,N1] = plane_orientation_from_circular_contours(data.K,data.T_cam,...
		data.isocontour.CurveParameters{1});
	[Xc2,N2] = plane_orientation_from_circular_contours(data.K,data.T_cam,...
		data.isocontour.CurveParameters{2});
	[L_1] = source_line_from_plane_pose(data.K,Xc1,psp1(1:3),psp1(4));
	[L_2] = source_line_from_plane_pose(data.K,Xc2,psp2(1:3),psp2(4));
	[ps,err] = closest_point_lines({L_1,L_2});
	dssp1 = distance_source_plane(psp1(1:3),ps,psp1(4));
	dssp2 = distance_source_plane(psp2(1:3),ps,psp2(4));
end
