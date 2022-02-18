function [dssp,psp,ps] = case_8_np_scaled(data)

	% Known parameters

	% We use data.groundtruth.ScenePlanePosition for scaling the scene, 
	%	so that our estimated plane pose psp{:} is close to those expected values

	% Estimation			
	i_p = 1;
	n_p = length(data.groundtruth.ScenePlaneDistanceSource);
	psp = cell(1,n_p);
	dssp = cell(1,n_p);
	centers = cell(1,n_p);
	all_planes = cell(2,n_p);
	all_im_center = cell(2,n_p);
	flag_ambig = 0;
	for i_p = 1:n_p
		pt_vis = get_visible_point_from_data(data,i_p);
		[Xc,N] = plane_orientation_from_circular_contours(data.K,data.T_cam,...
			data.isocontour.CurveParameters{i_p},pt_vis);
		if size(N,2)==2
			P_s_1 = source_plane_from_circle_center_orientation(data.K,Xc(:,1),N(:,1));
			P_s_2 = source_plane_from_circle_center_orientation(data.K,Xc(:,2),N(:,2));
			dssp{i_p} = {0};
			psp{i_p} = {[N(:,1);0],[N(:,2);0]};
			all_planes{1,i_p} = P_s_1;
			all_planes{2,i_p} = P_s_2;
			flag_ambig = 1;
		else
			P_s = source_plane_from_circle_center_orientation(data.K,Xc,N);
			dssp{i_p} = {0};
			psp{i_p} = {[N;0]};
			all_planes{1,i_p} = P_s;
			all_planes{2,i_p} = P_s;
		end
	end
	% If we have enough planes, we can calculate all their intersection as a line
	if n_p>1
		all_poss = 0:(2^(n_p)-1);
                err_poss = zeros(1,length(all_poss));
		L_poss = cell(1,length(all_poss));
		for i = 1:length(all_poss)
                        decomp = dec2bin(all_poss(i),n_p);
                        add_planes = cell(1,n_p);
                        for i_p=1:n_p
                                add_planes{i_p} = all_planes{str2num(decomp(i_p))+1,i_p};
                        end
			[L,err] = source_line_from_planes(add_planes);
                        L_poss{i} = L;
                        err_poss(i) = err;
                end
		% If we have more than 1 isocontour we have no ambiguities:
		% 	* The intersection of lines can be calculated
		%	* The scale can be adjusted based on the groundtruth
		if length(data.isocontour.CurveParameters{i_p})>1
			% There is no selection required (there are all the same values)
			% the minmum function selects the first one
			[~,i_min] = min(err_poss);
			ps = {L_poss{i_min}};
			% First suppose that the source is at a distance 1 of the camera
			[d,m] = plucker_matrix_to_dm(ps{1});
                        d = d/norm(d);
                        d = d*sign(d(3));
                        ps{1} = d;
                        % Then every other parameters can be calculated
                        for i_p = 1:n_p
                                % We calculate Xc again
                                [Xc,N] = plane_orientation_from_circular_contours(data.K,data.T_cam,...
                                        data.isocontour.CurveParameters{i_p});
                                dssp{i_p}{1} = distance_source_plane_from_circle_center_source(data.K,Xc,N,ps{1});
                                psp{i_p}{1}(4) = distance_plane_camera(N,ps{1},dssp{i_p}{1});
                        end
			% We can evaluate the distance ratio for each scene plane with groundtruth
			for i_p = 1:n_p
				data.groundtruth.ScenePlanePosition{i_p}
				psp{i_p}{1}(4) = 	
			end
		else
			% Otherwise we keep all the possible lines containing the source
			ps = L_poss;
		end
		if length(data.isocontour.CurveParameters{i_p})>1
			[d,m] = plucker_matrix_to_dm(ps{1});
			d = d/norm(d);
			d = d*sign(d(3));
			ps{1} = d*scale;
			% Then once the source is known, every other parameters can be calculated
			for i_p = 1:n_p
				% We calculate Xc again
				[Xc,N] = plane_orientation_from_circular_contours(data.K,data.T_cam,...
					data.isocontour.CurveParameters{i_p});
				dssp{i_p}{1} = distance_source_plane_from_circle_center_source(data.K,Xc,N,ps{1});
				psp{i_p}{1}(4) = distance_plane_camera(N,ps{1},dssp{i_p}{1});
			end
		end
	else
		if flag_ambig
			ps = transpose(all_planes{:,i_p});
		else
			ps = {all_planes{1,i_p}};
		end
	end
end
