function [dssp,psp,ps] = case_6_np_psp(data)

        % Known parameters
        n_p = length(data.groundtruth.ScenePlaneDistanceSource);
        i_p = 1;
        psp = cell(1,n_p);
	for i_p = 1:n_p
		psp{i_p} = {[data.groundtruth.ScenePlaneOrientation{i_p}(:,3);data.groundtruth.ScenePlanePosition{i_p}]};
	end

	% Estimation
	dssp = cell(1,n_p);
	all_lines = cell(2,n_p);
	flag_ambig = 0;
	for i_p = 1:n_p
		[Xc,N] = plane_orientation_from_circular_contours(data.K,data.T_cam,...
			data.isocontour.CurveParameters{i_p});
		if size(N,2)==2
			all_lines{1,i_p} = source_line_from_plane_pose(data.K,Xc(:,1),psp{i_p}{1}(1:3),psp{i_p}{1}(4));
			all_lines{2,i_p} = source_line_from_plane_pose(data.K,Xc(:,2),psp{i_p}{1}(1:3),psp{i_p}{1}(4));
			dssp{i_p} = {0};
			flag_ambig = 1;
		else
			all_lines{1,i_p} = source_line_from_plane_pose(data.K,Xc,psp{i_p}{1}(1:3),psp{i_p}{1}(4));
			all_lines{2,i_p} = source_line_from_plane_pose(data.K,Xc,psp{i_p}{1}(1:3),psp{i_p}{1}(4));
			dssp{i_p} = {0};
		end
	end
	if n_p>1
		all_poss = 0:(2^(n_p)-1);
                err_poss = zeros(1,length(all_poss));
                S_poss = zeros(3,length(all_poss));
                for i = 1:length(all_poss)
                        decomp = dec2bin(all_poss(i),n_p);
                        add_lines = cell(1,n_p);
                        for i_p=1:n_p
                                add_lines{i_p} = all_lines{str2num(decomp(i_p))+1,i_p};
                        end
                        [ps,err] = closest_point_lines(add_lines);
                        S_poss(:,i) = ps;
                        err_poss(:) = err;
                end
                [~,i_min] = min(err_poss);
                ps = {S_poss(:,i_min)};
                % Once we calculated the source we can actually calculate the distance now
                for i_p = 1:n_p
			dssp{i_p} = {distance_source_plane(psp{i_p}{1}(1:3),ps{1},psp{i_p}{1}(4))};
                end
	else
		ps = transpose(all_lines{:,1});
	end
end
