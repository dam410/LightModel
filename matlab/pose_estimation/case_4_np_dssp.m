function [dssp,psp,ps] = case_4_np_dssp(data)


	% Known parameters
        n_p = length(data.groundtruth.ScenePlaneDistanceSource);
        dssp = cell(1,n_p);
        for i_p = 1:n_p
                dssp{i_p} = {data.groundtruth.ScenePlaneDistanceSource{i_p}};
        end


	% Estimation
	i_p = 1;
	psp = cell(1,n_p);
	all_lines = cell(2,n_p);
	flag_ambig = 0;
	for i_p = 1:n_p
		[Xc,N] = plane_orientation_from_circular_contours(data.K,data.T_cam,...
			data.isocontour.CurveParameters{i_p});
		if size(N,2)==2
			% Much more complicated case, best solution would probably be to calculate
			%	the 2^{n_p} possible intersection point and keeping the best
			all_lines{1,i_p} = source_line_from_distance_source_plane(data.K,...
				Xc(:,1),N(:,1),dssp{i_p}{1});
			all_lines{2,i_p} = source_line_from_distance_source_plane(data.K,...
				Xc(:,2),N(:,2),dssp{i_p}{1});
			psp{i_p} = {[N(:,1);0],[N(:,2);0]};
			flag_ambig = flag_ambig+1;
		else
			all_lines{1,i_p} = source_line_from_distance_source_plane(data.K,Xc,N,dssp{i_p}{1});
			all_lines{2,i_p} = source_line_from_distance_source_plane(data.K,Xc,N,dssp{i_p}{1});
			psp{i_p} = {[N;0]};
		end
	end
	if n_p>1
		% Calculate the intersection of every possible n_p-uple of lines among 2^{n_p} possibilities
		all_poss = 0:(2^(n_p)-1);
		err_poss = zeros(1,length(all_poss));
		S_poss = zeros(3,length(all_poss));
		ps_temp = cell(1,length(all_poss));
		for i = 1:length(all_poss)
			decomp = dec2bin(all_poss(i),n_p);
			add_lines = cell(1,n_p);
			for i_p=1:n_p
				i_ambig = str2num(decomp(i_p))+1;
				add_lines{i_p} = all_lines{i_ambig,i_p};
			end
			[ps,err] = closest_point_lines(add_lines);
			S_poss(:,i) = ps;
			err_poss(i) = err;
			ps_temp{i} = ps;
		end
		if flag_ambig < n_p
			% Sometimes the best position is behind both plane we can remove source
			% that make the point behind it
			[~,i_min] = min(err_poss);
			ps = {S_poss(:,i_min)};
			% Once we calculated the source we can actually calculate the distance now
			decomp = dec2bin(all_poss(i_min),n_p);
			for i_p = 1:n_p
				i_ambig = str2num(decomp(i_p))+1;
				N = psp{i_p}{i_ambig}(1:3);
				d = distance_plane_camera(N,ps{1},dssp{i_p}{1});
				psp{i_p} = {[N;d]};
			end
		else
			ps = ps_temp;	
		end
	else
		ps = transpose(all_lines(:,1));
	end
end
