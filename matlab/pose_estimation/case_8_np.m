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
		% Get visible points if available
		pt_vis = get_visible_point_from_data(data,i_p);
		% Check if we have only one family of curves otherwise
		%	calculate Xc and N separately for both solutions in the same format
		%	and store the results	
		if iscell(data.isocontour.CurveParameters{i_p}{1})
			Xc = [];
			N = [];
			for i_ambig = 1:length(data.isocontour.CurveParameters{i_p})
				[Xc_ambig,N_ambig] = plane_orientation_from_circular_contours(data.K,data.T_cam,...
					data.isocontour.CurveParameters{i_p}{i_ambig},pt_vis);
				Xc = [Xc,Xc_ambig];
				N = [N,N_ambig];
			end
		else
			[Xc,N] = plane_orientation_from_circular_contours(data.K,data.T_cam,...
				data.isocontour.CurveParameters{i_p},pt_vis);
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
		all_poss = 0;
		err_poss = 0;
		if n_p>3
			[L,err] = source_line_from_planes_robust(all_planes(1,:));
		else
			[L,err] = source_line_from_planes(all_planes(1,:));
		end
		L_poss = {L};
		ps = L_poss;
		% If given a scale, we can calculate everything using the given scale
		if nargin>1 && length(data.isocontour.CurveParameters)>1
			% For each possible direction of the source, using the scale given by the source
			%	position, we can fix all the other parameters
			dssp = cell(1,n_p);

			% Calculate the position of the source
			[d,m] = plucker_matrix_to_dm(ps{1});
			d = d/norm(d);
			ps{1} = d*scale;
			% Check whether the point is in front or behidn the camera
			if isfield(data,'groundtruth')
				if norm(ps{1}+data.groundtruth.SourcePosition) < norm(ps{1}-data.groundtruth.SourcePosition)
				ps{1} = -ps{1};
			end
			% A solution is possible if any point X on the plane follows dot(N,S-X)>0
			for i_p = 1:n_p
				% Get a visible point
				pt_vis = get_visible_point_from_data(data,i_p);
				test_pos_shading = logical(zeros(1,length(psp{i_p})));
				for i_ambig = 1:length(psp{i_p})
					% Check shading equation on one visible point
					% Otherwise consider that brighest point is visible.
					% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
					% ! Warning if pt_vis is not given, the above is a strong assumption !
					% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
					% 	Recall the normal and center projection
					Xc = Xcs{i_p}(:,i_ambig);
					N = Ns{i_p}(:,i_ambig);
					% 	Calculate the relative distances with point source position
					[h,d] = distance_source_plane_from_circle_center_source_lines(...
						data.K,Xc,N,ps{1});
					% Check denegerate cases (scene plane is passing by cam center)
					if abs(d) < EPS
						h = Inf;
					end
					dssp{i_p}{i_ambig} = h;
					psp{i_p}{i_ambig} = [N;d];
					if ~isempty(pt_vis)
						X = [transpose(pt_vis);1];
						% Project X on the actual scene plane
						lambda_X = -d/dot(X,N);
						X = lambda_X*X;
					else
						X = ps{1} + dssp{i_p}{i_ambig}*psp{i_p}{i_ambig}(1:3);
					end
					% Check the shading equation, we use N pointing opposite way to cam
					test_pos_shading(i_ambig) = (dot(ps{1}-X,N)<0);
				end
				dssp{i_p} = dssp{i_p}(test_pos_shading);
				psp{i_p} = psp{i_p}(test_pos_shading);
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
