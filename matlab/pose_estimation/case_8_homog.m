function [dssp,psp,ps] = case_8_homog(data,scale)

	% Known parameters

	% Estimation			
	i_p = 1;
	n_p = length(data.groundtruth.ScenePlaneDistanceSource);
	psp = cell(1,n_p);
	dssp = cell(1,n_p);
	all_planes = cell(2,n_p);
	all_im_center = cell(2,n_p);
	flag_ambig = 0;
	% Saving the image of center and normals
	Xcs = cell(1,n_p);
	Ns = cell(1,n_p);
	for i_p = 1:n_p
		% Calculate the normal and the center using the homography
		nb_ambig = length(data.homography{i_p});
		N = zeros(3,nb_ambig);
		Xc = zeros(2,nb_ambig);
		for i_ambig = 1:nb_ambig
			H = data.homography{i_p}{i_ambig};
			invH = inv(H);
			[R,t] = decompose_homography(inv(data.K),H);
			N(:,i_ambig) = R(:,3);
			Xc_3D = invH*[0;0;1];
			Xc(:,i_ambig) = Xc_3D([1,2])/Xc_3D(3);
		end
		Xcs{i_p} = Xc;
		Ns{i_p} = N;
		if size(N,2)==2
			P_s_1 = source_plane_from_circle_center_orientation(data.K,Xc(:,1),N(:,1));
			P_s_2 = source_plane_from_circle_center_orientation(data.K,Xc(:,2),N(:,2));
			dssp{i_p} = {0,0};
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
	% If we have enough planes, we can calculate their intersection as a line
	if n_p>1
		if flag_ambig
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
			% If we have more than 2 planes we can solve the ambiguities
			%       the intersection of lines can be calculated
			if n_p>2 || ~flag_ambig
				[~,i_min] = min(err_poss);
				ps = {L_poss{i_min}};
				all_poss = all_poss(i_min);
			else
				all_poss = 0;
				ps = L_poss;
			end
		else
			[L,err] = source_line_from_planes(all_planes(1,:));
			ps = {L};
			all_poss = 0;
			err_poss = err;
		end
		% If given a scale, we calculate everything using the given scale
		if nargin>1
			nb_sol = length(all_poss);
			% Calculate the position of the source with the scale
			for i_ambig = 1:nb_sol
				[d,m] = plucker_matrix_to_dm(ps{i_ambig});
				d = d/norm(d);
				d = d*sign(d(3));
				ps{i_ambig} = d*scale;
			end
			% For each plane calculate all the possible poses
			for i_p = 1:n_p
				dssp{i_p} = cell(1,nb_sol);
				psp{i_p} = cell(1,nb_sol);
				for i_ambig = 1:nb_sol
					decomp_min = dec2bin(all_poss(i_ambig),n_p);
					desambig_indice = str2num(decomp_min(i_p))+1;
					% Once the source is known, every other parameters can be calculated
						% We calculate Xc again
					N = Ns{i_p}(:,desambig_indice);
					Xc = Xcs{i_p}(:,desambig_indice);
					dssp{i_p}{i_ambig} = distance_source_plane_from_circle_center_source(...
						data.K,Xc,N,ps{i_ambig});
					psp{i_p}{i_ambig} = [Ns{i_p}(:,desambig_indice);...
						distance_plane_camera(Ns{i_p}(:,desambig_indice),...
							ps{i_ambig},dssp{i_p}{i_ambig})];
				end
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
