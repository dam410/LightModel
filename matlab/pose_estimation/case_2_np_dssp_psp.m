function [dssp,psp,ps] = case_2_np_dssp_psp(data)

	% Known parameters
        n_p = length(data.groundtruth.ScenePlaneDistanceSource);
        dssp = cell(1,n_p);
        psp = cell(1,n_p);
        for i_p = 1:n_p
                dssp{i_p} = {data.groundtruth.ScenePlaneDistanceSource{i_p}};
                psp{i_p} = {[data.groundtruth.ScenePlaneOrientation{i_p}(:,3);data.groundtruth.ScenePlanePosition{i_p}]};
        end

	% Estimation 
	S_array = [];
	flag_ambig = 0;
	for i_p = 1:n_p
		[Xc,N] = plane_orientation_from_circular_contours(data.K,data.T_cam,...
			data.isocontour.CurveParameters{i_p});
		if size(N,2)==2
			S_1 = source_from_plane_pose_distance_source_plane(data.K,Xc(:,1),psp{i_p}{1}(1:3),psp{i_p}{1}(4),dssp{i_p}{1});
			S_2 = source_from_plane_pose_distance_source_plane(data.K,Xc(:,2),psp{i_p}{1}(1:3),psp{i_p}{1}(4),dssp{i_p}{1});
			S_array = [S_array;[S_1,S_2]];
			flag_ambig = 1;
		else
			S_1 = source_from_plane_pose_distance_source_plane(data.K,Xc,psp{i_p}{1}(1:3),psp{i_p}{1}(4),dssp{i_p}{1});
			S_array = [S_array;[S_1,S_1]];
		end
	end
	% Solving ambiguity using consensus
	if flag_ambig
		comp_vect =  @(x,y) norm(x-y);
		[V_best,index_best] = get_best_ambiguous(S_array,n_p,3,comp_vect);
		ps = {V_best};
	else
		ps = {mean(reshape(S_array(:,1),3,n_p),2)};
	end
end
