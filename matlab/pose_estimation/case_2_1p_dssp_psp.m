function [dssp,psp,ps] = case_2_1p_dssp_psp(data,i_surface)

	% Select the surface in the data 
	if nargin < 2
		i_surface = 1;
	end

	% Known parameters
	dssp = {{data.groundtruth.ScenePlaneDistanceSource{i_surface}}};
	psp = {{[data.groundtruth.ScenePlaneOrientation{i_surface}(:,3);...
		data.groundtruth.ScenePlanePosition{i_surface}]}};
	i_p = 1;
	% Estimation 
	[Xc,N] = plane_orientation_from_circular_contours(data.K,data.T_cam,...
		data.isocontour.CurveParameters{i_surface});
	if size(N,2)==2
		ps = {...
			source_from_plane_pose_distance_source_plane(data.K,...
				Xc(:,1),psp{i_p}{1}(1:3),psp{i_p}{1}(4),dssp{i_p}{1}),...
			source_from_plane_pose_distance_source_plane(data.K,...
				Xc(:,2),psp{i_p}{1}(1:3),psp{i_p}{1}(4),dssp{i_p}{1}),...
		};
	else
		ps = {source_from_plane_pose_distance_source_plane(data.K,...
			Xc,psp{i_p}{1}(1:3),psp{i_p}{1}(4),dssp{i_p}{1})};
	end
end
