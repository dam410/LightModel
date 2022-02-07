% This function will be more generic and can be used for any number of planes and estimation
function [err_h,err_orient,err_d,err_s] = eval_pose_n(dssp,psp,ps,data,i_surface)
	

	n_plane = length(dssp);

	% Error for the distance estimation between scene plane and source
	err_h = zeros(1,n_plane);
	for i_p = 1:n_plane

		if nargin < 5
			i_surface = i_p;
		end

		% Calculate the minimum, among all the solutions
		%	If no solution is proposed dssp{i_p} error is infinite
		if ~isempty(cell2mat(dssp{i_p}))
			err_h(i_p) = min(abs(cell2mat(dssp{i_p})-data.groundtruth.ScenePlaneDistanceSource{i_surface}));
		else
			err_h(i_p) = Inf;
		end
	end

	% Error for the orientation of the scene plane
	err_orient = zeros(1,n_plane);
	for i_p = 1:n_plane

		if nargin < 5
			i_surface = i_p;
		end

		% Calculate the minimum among all the solutions
		n_sol = length(psp{i_p});
		Planes_sol = cell2mat(psp{i_p});
		%data.groundtruth.ScenePlaneOrientation{i_surface}(:,3)*ones(1,n_sol)
		if ~isempty(Planes_sol)
			scalar_prod = dot(data.groundtruth.ScenePlaneOrientation{i_surface}(:,3)*ones(1,n_sol),Planes_sol(1:3,:));
			err_orient(i_p) = min(180/pi*acos(min(1,max(scalar_prod,0))));
		else
			err_orient(i_p) = Inf;
		end
	end

	% Error for the distance between origin and scene plane
	err_d = zeros(1,n_plane);
	for i_p = 1:n_plane

		if nargin < 5
			i_surface = i_p;
		end
		n_sol = length(psp{i_p});
		Planes_sol = cell2mat(psp{i_p});
		if ~isempty(Planes_sol)
			err_d(i_p) = min(abs(Planes_sol(4,:)-data.groundtruth.ScenePlanePosition{i_surface}));
		else
			err_d(i_p) = Inf;
		end
	end
	
	% The evaluation of the source error depends on what we calculated
	if size(ps{1},1) == 3 && size(ps{1},2) == 1
		% We have a 3D point
		diff_s = cell2mat(ps)-data.groundtruth.SourcePosition;
		err_s = min(sqrt(sum(diff_s.^2,1)));
	elseif size(ps{1},1) == 4 && size(ps{1},2) == 1
		% We have a plane
		% Check normalization of the plane
		planes_s = cell2mat(ps);
		planes_s = planes_s./sqrt(sum(planes_s(1:3,:).^2,1));
		% Distance of the point to the plane
		err_s = min(abs(dot(planes_s(1:3,:),data.groundtruth.SourcePosition*ones(1,length(ps)))+planes_s(4,:)));
	elseif size(ps{1},1)==4 && size(ps{1},2) == 4
		if length(ps)>0
			err_s_array = zeros(1,length(ps));
			for i_s=1:length(ps)
				% We have a line in plucker matrix coordinates
				[d,m] = plucker_matrix_to_dm(ps{i_s});
				% Using the reciprocal product
				err_s_array(i_s) = norm(cross_antisym(d)*data.groundtruth.SourcePosition-m);
			end
			err_s = min(err_s_array);
		else
			err_s = Inf;
		end
	else
		disp('Not recognized source output');
		err_s = 0;
	end
end
