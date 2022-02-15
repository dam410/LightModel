function [dssp,psp,ps] = case_8_np(data,scale)

	EPS = 1e-8;

	% Known parameters

	% Estimation			
	i_p = 1;
	n_p = length(data.isocontour.CurveParameters);
	psp = cell(1,n_p);
	dssp = cell(1,n_p);
	centers = cell(1,n_p);
	all_planes = cell(2,n_p);
	all_im_center = cell(2,n_p);
	flag_ambig = 0;
	Xcs = cell(1,n_p);
	Ns = cell(1,n_p);
	for i_p = 1:n_p
		% Check if we have only one family of curves otherwise
		%	calculate Xc and N separately for both solutions in the same format
		%	and store the results	
		if iscell(data.isocontour.CurveParameters{i_p}{1})
			Xc = [];
			N = [];
			for i_ambig = 1:length(data.isocontour.CurveParameters{i_p})
				[Xc_ambig,N_ambig] = plane_orientation_from_circular_contours(data.K,data.T_cam,...
					data.isocontour.CurveParameters{i_p}{i_ambig});
				Xc = [Xc,Xc_ambig];
				N = [N,N_ambig];
			end
		else
			[Xc,N] = plane_orientation_from_circular_contours(data.K,data.T_cam,...
				data.isocontour.CurveParameters{i_p});
		end
		Xcs{i_p} = Xc;
		Ns{i_p} = N;
		if size(N,2)==2
			P_s_1 = source_plane_from_circle_center_orientation(data.K,Xc(:,1),N(:,1));
			P_s_2 = source_plane_from_circle_center_orientation(data.K,Xc(:,2),N(:,2));
			%disp('N Xc Ps angle');
			%N(:,1)
			%Xc(:,1)
			%P_s_1(:,1)
			%180*acos([1,0,0]*P_s_1(1:3,1))/pi
			dssp{i_p} = {Inf};
			psp{i_p} = {[N(:,1);0],[N(:,2);0]};
			all_planes{1,i_p} = P_s_1;
			all_planes{2,i_p} = P_s_2;
			flag_ambig = 1;
		else
			P_s = source_plane_from_circle_center_orientation(data.K,Xc,N);
			dssp{i_p} = {Inf};
			psp{i_p} = {[N;0]};
			all_planes{1,i_p} = P_s;
			all_planes{2,i_p} = P_s;
		end
	end
	% If we have enough planes, we can calculate their intersection as a line
	if n_p>1
		% After some check indeed the vectors in each pairs will be opposite anyway so no need to compute
		%	all the posssibility just consider one case to calculate ps
		%if flag_ambig
		%	all_poss = 0:(2^(n_p)-1);
		%	err_poss = zeros(1,length(all_poss));
		%	L_poss = cell(1,length(all_poss));
		%	for i = 1:length(all_poss)
		%		decomp = dec2bin(all_poss(i),n_p);
		%		add_planes = cell(1,n_p);
		%		for i_p=1:n_p
		%			add_planes{i_p} = all_planes{str2num(decomp(i_p))+1,i_p};
		%		end
		%		[L,err] = source_line_from_planes(add_planes);
		%		L_poss{i} = L;
		%		err_poss(i) = err;
		%	end
		%	display_normals(all_planes);
		%else
			%display_normals(all_planes);
		all_poss = 0;
		err_poss = 0;
		if n_p>3
			[L,err] = source_line_from_planes_robust(all_planes(1,:));
		else
			[L,err] = source_line_from_planes(all_planes(1,:));
		end
		L_poss = {L};
		%end
		% Does not solve the ambiguity with intersection of the lines
		%	
		ps = L_poss;
		%all_planes{1,1}
		%all_planes{2,1}
		%all_planes{1,2}
		%all_planes{2,2}
		%all_planes{1,4}
		%% If we have more than 1 isocontour we have no ambiguities
		%% 	the intersection of lines can be calculated
		%if length(data.isocontour.CurveParameters{i_p})>1
		%	[~,i_min] = min(err_poss);
		%	ps = {L_poss{i_min}};
		%else
		%	ps = L_poss;
		%end
		% If given a scale, we can calculate everything using the given scale
		if nargin>1 && length(data.isocontour.CurveParameters)>1
			% For each possible direction of the source, using the scale given by the source
			%	position, we can fix all the other parameters
			dssp = cell(1,n_p);

			% Calculate the position of the source
			[d,m] = plucker_matrix_to_dm(ps{1});
			d = d/norm(d);
			d = d*sign(d(3));
			ps{1} = d*scale;

			% This is not possible anymore as the source is determined but this prevent from knowing the ambiguity for the plane
			%% Then once the source is known, every other parameters can be calculated
			%for i_p = 1:n_p
			%	dssp{i_p} = cell(1,length(ps));
			%	psp{i_p} = cell(1,length(ps));
			%	for i_ambig = 1:length(ps)
			%		decomp_min = dec2bin(all_poss(i_ambig),n_p);
			%		desambig_indice = str2num(decomp_min(i_p))+1;
			%		if size(ps{i_ambig}) == [4,4]
			%			[d,m] = plucker_matrix_to_dm(ps{i_ambig});
			%			d = d/norm(d);
			%			d = d*sign(d(3));
			%			ps{i_ambig} = d*scale;
			%		end
			%		% Recall the normal and center projection
			%		Xc = Xcs{i_p}(:,desambig_indice);
			%		N = Ns{i_p}(:,desambig_indice);
			%		% Instead of using this function we may use in fact the intersection between
			%		% normal line to the scene plane and back projection of circle center
			%		%dssp{i_p}{1} = distance_source_plane_from_circle_center_source(data.K,Xc,N,ps{1});
			%		% Cool but useless does not increase results
			%		[h,d] = distance_source_plane_from_circle_center_source_lines(...
			%			data.K,Xc,N,ps{i_ambig});
			%		dssp{i_p}{i_ambig} = h;
			%		psp{i_p}{i_ambig} = [N;d];
			%	end
			%end

			% Now we can check if there is no case where the visible point is placed behind the visible
			%	part of the plane.
			% Cond 1 : If a point X is visible by the camera with a positive shading : N^{T}(O-X)>0
			% 	then all points of the plane which are visible by the camera will have a positive
			%	shading, so we only have to check on one specific point, the hightlight point.
			% Only do that if there are more than one solution
			if length(ps)>1

				% Note that this case should not happen anymore
				test_pos_shading = ones(1,length(ps));
				for i_ambig = 1:length(ps)
					for i_p = 1:n_p
						% Calculate the position of the point
						X = ps{i_ambig} + dssp{i_p}{i_ambig}*psp{i_p}{i_ambig}(1:3);
						test_pos_shading(i_ambig) = (dot(ps{i_ambig}-X,-X)>0)*test_pos_shading(i_ambig);
					end
				end
				ps = ps(find(test_pos_shading));
				for i_p = 1:n_p
					dssp{i_p} = dssp{i_p}(find(test_pos_shading));
					psp{i_p} = psp{i_p}(find(test_pos_shading));
				end
			% Only one source but for each ambiguity we check if its plane can be physically lit
			%	and calculate its full pose as well
			else
				for i_p = 1:n_p
					test_pos_shading = logical(zeros(1,length(psp{i_p})));
					for i_ambig = 1:length(psp{i_p})
						% Recall the normal and center projection
						Xc = Xcs{i_p}(:,i_ambig);
						N = Ns{i_p}(:,i_ambig);
						% Instead of using this function we may use in fact the intersection between
						% normal line to the scene plane and back projection of circle center
						%dssp{i_p}{1} = distance_source_plane_from_circle_center_source(data.K,Xc,N,ps{1});
						% Cool but useless does not increase results
						[h,d] = distance_source_plane_from_circle_center_source_lines(...
							data.K,Xc,N,ps{1});
						if abs(d) < EPS
							h = Inf;
						end
						dssp{i_p}{i_ambig} = h;
						psp{i_p}{i_ambig} = [N;d];
						X = ps{1} + dssp{i_p}{i_ambig}*psp{i_p}{i_ambig}(1:3);
						test_pos_shading(i_ambig) = (dot(ps{1}-X,-X)>0);
						% If we have access to the points, we can check that 
						%	their reprojection considering plane pose are in front of the cam
						if isfield(data.isocontour,'Points')
							% Calculate for one isocontour point its real position in space
							%	Keep old format compatibility (without ambiguity in points)
							if iscell(data.isocontour.Points{i_p}{1})
								X_iso = [transpose(data.isocontour.Points{i_p}{1}{1}(1,:));1];
							else
								X_iso = [transpose(data.isocontour.Points{i_p}{1}(1,:));1];
							end
							X_iso = inv(data.K)*X_iso;
							lambda = -psp{i_p}{i_ambig}(4)/(transpose(psp{i_p}{i_ambig}(1:3))*X_iso);
							% Keep the solution if it gives points in front of the camera
							% 	But also if the the source is colocalised with the camera
							test_pos_shading(i_ambig) = (test_pos_shading(i_ambig) & lambda>0) | (norm(ps{1}) < EPS);
						end
					end
					dssp{i_p} = dssp{i_p}(test_pos_shading);
					psp{i_p} = psp{i_p}(test_pos_shading);
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
