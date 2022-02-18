function [dssp,psp,ps] = case_7_np_ps(data)

	% Known parameters
	ps = {data.groundtruth.SourcePosition};
	K = data.K;
	% Estimation
	if isfield(data,'isocontour')
		n_p = length(data.isocontour.CurveParameters);
	elseif isfield(data,'homography')
		n_p = length(data.isocontour.homography);
	else
		n_p = 1;
	end
	dssp = cell(1,n_p);
	psp = cell(1,n_p);
	for i_p = 1:n_p
		pt_vis = get_visible_point_from_data(data,i_p);
		% Check if we have only one family of curves otherwise
                %       calculate Xc and N separately for both solutions in the same format
                %       and store the results
                if iscell(data.isocontour.CurveParameters{i_p}{1})
                        Xc = [];
                        N = [];
                        for i_ambig = 1:length(data.isocontour.CurveParameters{i_p})
                                [Xc_ambig,N_ambig] = plane_orientation_from_circular_contours(data.K,data.T_cam,...
                                        data.isocontour.CurveParameters{i_p}{i_ambig},pt_vis);
                                Xc = [Xc,Xc_ambig];
                                N = [N,N_ambig];
                        end
			N = -N;

                else
                        [Xc,N] = plane_orientation_from_circular_contours(data.K,data.T_cam,...
                                data.isocontour.CurveParameters{i_p},pt_vis);
                end
		if size(N,2)==2
			% NO IT IS NOT WORKING, THE AMBIGUOUS POSE OF THE PLANE GIVE THE SAME PLANE !!!!!!!!!!!!!
			%% The ambiguity can be solved in this case
			%% 	The normal and the ray back projection of the circle center formed a plane
			%%	This plane passing through the camera center must also intersect the source
			%%	So a criteria would be to select the pair (Xc,N) so that the image of the source
			%%	by the camera is clostest to the image of this plane ( a line )
			%P_1 = source_plane_from_circle_center_orientation(data.K,Xc(:,1),N(:,1));
			%P_2 = source_plane_from_circle_center_orientation(data.K,Xc(:,2),N(:,2));
			%dist_planes = [dot(P_1(1:3),ps{1})+P_1(4),dot(P_2(1:3),ps{1})+P_2(4)]
			%[~,i_ambig] = min(abs(dist_planes));
			%dssp{i_p} = {distance_source_plane_from_circle_center_source(data.K,Xc(:,i_ambig),N(:,i_ambig),ps{1})};
			%d = distance_plane_camera(N(:,i_ambig),ps{1},dssp{i_p}{1});
			%psp{i_p} = {[N(:,i_ambig);d]};
			h = [distance_source_plane_from_circle_center_source(K,Xc(:,1),N(:,1),ps{1}),...
				distance_source_plane_from_circle_center_source(K,Xc(:,2),N(:,2),ps{1})];
			d = [distance_plane_camera(N(:,1),ps{1},h(1)),distance_plane_camera(N(:,2),ps{1},h(2))];
			dir_S = [dir_S_X(K,ps{1},Xc(:,1),N(:,1),d(1)),dir_S_X(K,ps{1},Xc(:,2),N(:,2),d(2))];
			% TODO Do not know why it is working here but it works with a minus
			[~,i_ambig] = max(-dot(dir_S,N));
			psp{i_p} = {[N(:,i_ambig);d(i_ambig)]};
			dssp{i_p} = {h(i_ambig)};
		else
			dssp{i_p} = {distance_source_plane_from_circle_center_source(K,Xc,N,ps{1})};
			d = distance_plane_camera(N,ps{1},dssp{i_p}{1})
			psp{i_p} = {[N;d]};
		end
	end
end
