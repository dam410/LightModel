function [dssp,psp,ps] = case_3_1p_dssp_ps(data)

	% Known parameters
	K = data.K;%*[0,1,0;1,0,0;0,0,1];
	n_p = length(data.groundtruth.ScenePlaneDistanceSource);
	dssp = cell(1,n_p);
	psp = cell(1,n_p);
	for i_p = 1:n_p
		dssp{i_p} = {data.groundtruth.ScenePlaneDistanceSource{i_p}};
	end
	ps = {data.groundtruth.SourcePosition};
	
	% Estimation
	psp = cell(1,n_p);
	for i_p = 1:n_p
		[Xc,N] = plane_orientation_from_circular_contours(K,data.T_cam,...
			data.isocontour.CurveParameters{i_p});
		if size(N,2)==2
			d_1 = distance_plane_camera(N(:,1),ps{1},dssp{i_p}{1});
			d_2 = distance_plane_camera(N(:,2),ps{1},dssp{i_p}{1});
			psp{i_p} = {[N(:,1);d_1],[N(:,2);d_2]};
		else
			d = distance_plane_camera(N,ps{1},dssp{i_p}{1});
			psp{i_p} = {[N;d]};
		end
	end
end
