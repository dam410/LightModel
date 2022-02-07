function [dssp,psp,ps] = case_3_1p_dssp_ps(data,i_surface)

	% Select the surface in the data 
        if nargin < 2
                i_surface = 1;
        end

	% Known parameters
	dssp = {{data.groundtruth.ScenePlaneDistanceSource{i_surface}}};
	ps = {data.groundtruth.SourcePosition};
	i_p = 1;
	% Estimation
	[Xc,N] = plane_orientation_from_circular_contours(data.K,data.T_cam,...
		data.isocontour.CurveParameters{i_surface});
	if size(N,2)==2
		d_1 = distance_plane_camera(N(:,1),ps{1},dssp{i_p}{1});
		d_2 = distance_plane_camera(N(:,2),ps{1},dssp{i_p}{1});
		psp = {{[N(:,1);d_1],[N(:,2);d_2]}};
	else
		d = distance_plane_camera(N,ps{1},dssp{i_p}{1});
		psp = {{[N;d]}};
	end
end
