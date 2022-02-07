function [dssp,psp,ps,err_opt] = refine_pose_isocontours(data,dssp,psp,ps,speed,constraints)
% We will suppose first that everyhting has been estimated function 
	% Number of planes
	n_p = length(data.isocontour.CurveParameters);
	radius = cell(1,n_p);
	ang_pts = cell(1,n_p);
	n_iso = length(data.isocontour.CurveParameters{1});
	cell_arg = cell(1,5+2*n_p*n_iso);
	cell_arg{1} = data.K;
	cell_arg{2} = ps{1};
	cell_arg{3} = zeros(3,n_p);
	cell_arg{4} = zeros(1,n_p);
	cell_arg{5} = zeros(n_iso,n_p);
	i_cur = 6;
	data_copy = data;
	pts_2D_data = [];
	NB_MAX_PTS = 100;
	if nargin < 6
		constraints = 'None';
	end
	% Initialize the solution
	for  i_p=1:n_p
		% Treat the case when no solution have been found
		if isempty(psp{i_p})
			X = ps{1} - [0;0;1];
			N = [0;0;1];
			h = -(-1+[0,0,1]*ps{1});
			d = -1;
		else
			X = ps{1} + dssp{i_p}{1}*psp{i_p}{1}(1:3);
			N = psp{i_p}{1}(1:3);
			h = dssp{i_p}{1};
			d = psp{i_p}{1}(4);
		end
		try
			R_plane = [null(N*transpose(N)),N];
		catch
			R_plane = eye(3);
		end
		if det(R_plane) < 0
			R_plane(:,1) = -R_plane(:,1);
		end
		R_plane = eul2rotm(rotm2eul(R_plane));
		cell_arg{3}(:,i_p) = transpose(rotm2eul(R_plane));
		cell_arg{4}(i_p) = h;
		if iscell(data_copy.isocontour.CurveParameters{i_p}{1})
			curveParameters = data_copy.isocontour.CurveParameters{i_p}{1};
			isoPoints = data_copy.isocontour.Points{i_p}{1};
		else
			curveParameters = data_copy.isocontour.CurveParameters{i_p};
			isoPoints = data_copy.isocontour.Points{i_p};
		end
		for i_iso = 1:n_iso
			ell = curveParameters{i_iso};
			data_copy.isocontour.CurveParameters{i_p}{i_iso} = ell;
			E = param2ellipse(ell);
			% Calculate the position of the estimated circle center
			% Calculate the back projection of the points on the plane
			% Always keep 500 points for each isocontours -> This means we can calculate the size of the vectors
			nb_pts_total = length(isoPoints{i_iso});
			if  nb_pts_total>NB_MAX_PTS
				samples = datasample(1:nb_pts_total,NB_MAX_PTS,'Replace',false);
				data_copy.isocontour.Points{i_p}{i_iso} = isoPoints{i_iso}(samples,:);
			else
				samples = datasample(1:nb_pts_total,NB_MAX_PTS,'Replace',true);
                                data_copy.isocontour.Points{i_p}{i_iso} = isoPoints{i_iso}(samples,:);
			end
			pts_2D_data = [pts_2D_data;data_copy.isocontour.Points{i_p}{i_iso}];
			nb_pts = length(data_copy.isocontour.Points{i_p}{i_iso});
			%plot(data_copy.isocontour.Points{i_p}{i_iso}(:,2),data_copy.isocontour.Points{i_p}{i_iso}(:,1),'+r');
			P_pts_ = [transpose(data_copy.isocontour.Points{i_p}{i_iso}(:,[2,1]));...
				ones(1,nb_pts)];
			P_pts_ = inv(data_copy.K)*P_pts_;
			lambda_pts = -d./(transpose(N)*P_pts_);
			P_pts = lambda_pts.*P_pts_;
			Y = transpose(R_plane)*X;
			P_plan = transpose(R_plane)*P_pts;
			P_plan_pts = P_plan-Y;
			r_pts = sqrt(P_plan_pts(1,:).^2+P_plan_pts(2,:).^2);
			ang_pts = atan2(P_plan_pts(2,:),P_plan_pts(1,:));
			cell_arg{5}(i_iso,i_p) = mean(r_pts);
			cell_arg{i_cur} = data_copy.isocontour.Points{i_p}{i_iso};
			cell_arg{i_cur+1} = ang_pts;
			i_cur = i_cur+2;

			%figure('Name','Display the circle and Y');
			%plot3(P_pts(1,:),P_pts(2,:),P_pts(3,:),'+r');
			%axis equal;
	
			%figure('Name','Display the circle and Y');
			%plot3(P_plan(1,:),P_plan(2,:),P_plan(3,:),'+r');
			%hold on;
			%plot3(Y(1),Y(2),Y(3),'+b');
			%axis equal;
		
			%figure('Name','Radius calculated');
			%plot(r_pts);

			%% Display the reprojection in 3D
			pts = isoPoints{i_iso};
			angles = ang_pts;
			pts_3D = X+R_plane*[mean(r_pts)*cos(angles);mean(r_pts)*sin(angles);zeros(1,length(angles))];
			pts_2D = data_copy.K*pts_3D;
			pts_2D = pts_2D./pts_2D(3,:);
			pts_detected = transpose(pts);

			%figure('Name','Display the circle and Y');
			%plot(pts_2D(1,:),pts_2D(2,:),'+r');
			%hold on;
			%plot(pts_detected(2,:),pts_detected(1,:),'+b');
		end
	end
	% Optimize the function now using our initial set of parameters
	%f_opt = @(parameters_vectors) function_to_optimized(K,parameters_vectors(1:
	[err] = function_to_optimized(cell_arg{:});
	[x_init,J_sparse] = cell_to_param_vect(cell_arg{:});
	%fun_opt = @(param) subfunction_proxy(param,data_copy.K,n_iso,n_p,data_copy);
	switch constraints
		case 'None'
			fun_opt_vect = @(param) vectorized_function_opt(data_copy.K,n_iso,n_p,pts_2D_data,NB_MAX_PTS,param);
		case 'Source'
			x_init_0 = x_init;
			fun_opt_vect = @(param) vectorized_function_opt(...
				data_copy.K,n_iso,n_p,pts_2D_data,NB_MAX_PTS,[x_init_0(1:3);param]);
			x_init = x_init(4:end);
		case 'PlanesOrientation'
			x_init_0 = x_init;
			fun_opt_vect = @(param) vectorized_function_opt(...
				data_copy.K,n_iso,n_p,pts_2D_data,NB_MAX_PTS,[param(1:3);x_init_0(4:(4+(3*n_p-1)));param(4:end)]);
			x_init(4:(4+(3*n_p-1))) = [];
		case 'PlanesPose'
			x_init_0 = x_init;
			h_vector = transpose(x_init_0((4+(3*n_p)):(3+4*n_p)));
			Rs = eul2rotm(transpose(reshape(x_init_0(4:(4+(3*n_p-1))),3,n_p)));
			d_vector = -h_vector-dot(x_init_0(1:3)*ones(1,n_p),reshape(Rs(:,3,:),3,n_p));
			fun_opt_vect = @(param) vectorized_function_opt_d_parameters(...
				data_copy.K,n_iso,n_p,pts_2D_data,NB_MAX_PTS,[param(1:3);x_init_0(4:(4+(3*n_p-1)));transpose(d_vector);param(4:end)]);
			% We have to change the h_vector into a d_vector
			x_init(4:(3+4*n_p)) = [];
		otherwise
			fun_opt_vect = @(param) vectorized_function_opt(data_copy.K,n_iso,n_p,pts_2D_data,NB_MAX_PTS,param);
	end
	err_vect = fun_opt_vect(x_init);
	options = optimoptions('lsqnonlin');
	if nargin>4 && speed==1
		options.Algorithm = 'trust-region-reflective';
		options.Display = 'off';
		options.UseParallel = true;
		options.JacobPattern = J_sparse;
		tic;
		[x_opt,err_opt] = lsqnonlin(fun_opt_vect,x_init,[],[],options);
		toc
	elseif nargin>4 && speed==2
		tic;
		fun_scalar = @(param) norm(fun_opt_vect(param));
		[x_opt,err_opt] = fminunc(fun_scalar,x_init);
                toc
	elseif nargin>4 && speed ==3
		% Debug function t understand what's wrong with outliers points
		fun_opt_vect = @(param) vectorized_function_opt(data_copy.K,n_iso,n_p,pts_2D_data,NB_MAX_PTS,param,'on');
		x_opt = x_init;
		err_opt = fun_opt_vect(x_init);
	else
		options.Algorithm = 'levenberg-marquardt';
		options.Display = 'off';
		options.UseParallel = true;
		options.JacobPattern = J_sparse;
		options.MaxIterations = 30;
                tic;
                [x_opt,err_opt] = lsqnonlin(fun_opt_vect,x_init,[],[],options);
                toc
	end
	switch constraints
		case 'None'
			%% Display for debugging
			%vectorized_function_opt(data_copy.K,n_iso,n_p,pts_2D_data,NB_MAX_PTS,x_init,'display_on');
			%vectorized_function_opt(data_copy.K,n_iso,n_p,pts_2D_data,NB_MAX_PTS,x_opt,'display_on');
			[dssp,psp,ps] = vector_to_pose(x_opt,data_copy.K,n_iso,n_p,data_copy);
		case 'Source'
			[dssp,psp,ps] = vector_to_pose([ps{1};x_opt],data_copy.K,n_iso,n_p,data_copy);
		case 'PlanesOrientaton'
			[dssp,psp,ps] = vector_to_pose([x_opt(1:3);x_init_0(4:(4+(3*n_p-1)));x_opt(4:end)],...
				data_copy.K,n_iso,n_p,data_copy);
		case 'PlanesPose'
			% Reclaculate h_vector from the d_vector
			h_vector = -d_vector-dot(x_opt(1:3)*ones(1,n_p),reshape(Rs(:,3,:),3,n_p));
			[dssp,psp,ps] = vector_to_pose([x_opt(1:3);x_init_0(4:(4+(3*n_p-1)));transpose(h_vector);x_opt(4:end)],...
				data_copy.K,n_iso,n_p,data_copy);
		otherwise
			[dssp,psp,ps] = vector_to_pose(x_opt,data_copy.K,n_iso,n_p,data_copy);
	end
	err_opt = fun_opt_vect(x_opt);
end

% Get back the parameters to original pose format
function [dssp,psp,ps] = vector_to_pose(vector,K,n_iso,n_p,data)
	cell_arg = param_vect_to_cell(vector,K,n_iso,n_p,data);

        angles_euler = cell_arg{3};
        h_vector = cell_arg{4};

	ps = cell(1,1);
	dssp = cell(1,n_p);
	psp = cell(1,n_p);
        ps{1} = cell_arg{2};

        for n=1:n_p
                R = eul2rotm(transpose(angles_euler(:,n)));
		d = distance_plane_camera(R(:,3),ps{1},h_vector(n));
		dssp{n} = {h_vector(n)};
		psp{n} = {[R(:,3);d]};
        end
end


function [err] = subfunction_proxy(vector,K,n_iso,n_p,data)
	%disp('Vectorization :');
	%tic;
	cell_arg = param_vect_to_cell(vector,K,n_iso,n_p,data);
	%time_passed = toc();
	%disp(['Temps : ',num2str(time_passed)]);	
	%disp('Fonction réelle :');
	%tic;
	err = function_to_optimized(cell_arg{:});
	%time_passed = toc();
	%disp(['Temps : ',num2str(time_passed)]);	
end

% Function must be call using function_to_optimized(S,angles_euler,d_vector,pts_1,angles_1,...pts_n,angles_n)
%	K : 3x3 matrix
%	S : 3x1 matrix
%	angles_euler : 3xn matrix
%	d_vector : 1xn matrix
%	r_vector : mxn matrix (Constraint the fact that there are the same number of isocontours for each plane)
%	pts_(1,1) : 2x(k_(1,1)) matrix
%	angles_(1,2) : 1x(k_(1,2)) matrix
%	...
%	pts_(1,m) : 2x(k_(1,m)) matrix
%	angles_(1,m) : 1x(k_(1,m)) matrix
%	...
%	pts_(n,1) : 2x(k_(n,1)) matrix
%	angles_(n,1) : 1x(k_(n_1)) matrix
%	...
%	pts_(n,m) : 2x(k_(m,n)) matrix
%	angles_(n,m) : 1x(k_(m,n)) matrix
function [err] = function_to_optimized(varargin)
	K = varargin{1};
	S = varargin{2};
	angles_euler = varargin{3};
	h_vector = varargin{4};
	n_plane = size(angles_euler,2);
	r_vector = varargin{5};
	m_iso = size(r_vector,1);
	err = [];
	for n=1:n_plane
		R = eul2rotm(transpose(angles_euler(:,n)));
		X = S + h_vector(n)*R(:,3);
		for m=1:m_iso
			i_current = 6+(n-1)*m_iso*2+(m-1)*2;
			pts = varargin{i_current};
			angles = varargin{i_current+1};
			pts_3D = X+R*[r_vector(m,n)*cos(angles);r_vector(m,n)*sin(angles);zeros(1,length(angles))];
			pts_2D = K*pts_3D;
			pts_2D = pts_2D./pts_2D(3,:);
			pts_2D_x = pts_2D(2,:);
			pts_2D_y = pts_2D(1,:);
			pts_x = transpose(pts(:,1));
			pts_y = transpose(pts(:,2));
			err_loc = [pts_x-pts_2D_x;pts_y-pts_2D_y];
			err = [err;err_loc(:)];
			%figure('Name','Display the reprojected points');
			%plot(pts_2D_x,pts_2D_y,'+r');
			%hold on;
			%plot(pts_x,pts_y,'+b');
			%figure('Name','Show the error');
			%plot(err_loc,'+r');
		end
	end
end

function [err] = vectorized_function_opt(K,n_iso,n_p,pts_2D,nb_pts,param_vect,display)
	S = param_vect(1:3);
	Rs = eul2rotm(transpose(reshape(param_vect(4:(4+(3*n_p-1))),3,n_p)));
	h_vector = transpose(param_vect((4+(3*n_p)):(3+4*n_p)));
	r_vector = reshape(param_vect((4+4*n_p):(3+(n_iso+4)*n_p)),n_iso,n_p);
	Xs = S+h_vector.*reshape(Rs(:,3,:),3,n_p);
	angles = reshape(param_vect((4+(n_iso+4)*n_p:end)),nb_pts*n_iso,n_p);
	cos_angles = cos(angles);
	sin_angles = sin(angles);
	zeros_angles = zeros(size(angles));
	err = zeros(n_iso*n_p*2*nb_pts,1);
	if nargin>6
		% Displaying the error
		figure('Name','Difference between isocontour points and their model prediction');
		imshow(ones(2072,4608));
		hold on;
	end
	for i_p=1:n_p
		% Calculate if the light source is on the same side
		%	of the plane to the optical center
		d = -h_vector(i_p)-transpose(S)*Rs(:,3,i_p);
		
		r_angles = repmat(transpose(r_vector(:,i_p)),nb_pts,1);
		X_pts = [r_angles(:).*cos_angles(:,i_p),r_angles(:).*sin_angles(:,i_p),zeros_angles(:,i_p)]*transpose(Rs(:,:,i_p));
		X_pts_2D = K*(Xs(:,i_p)+transpose(X_pts));
		X_pts_2D = X_pts_2D([2,1],:)./X_pts_2D(3,:);
		pts_2D_i_p = pts_2D((1+(i_p-1)*n_iso*nb_pts):i_p*n_iso*nb_pts,:);
		diff = pts_2D_i_p-transpose(X_pts_2D);
		if nargin>6
			plot(transpose([pts_2D_i_p(:,2),transpose(X_pts_2D(2,:))]),transpose([pts_2D_i_p(:,1),transpose(X_pts_2D(1,:))]),'-b');
			plot(pts_2D_i_p(:,2),pts_2D_i_p(:,1),'+g');
			plot(X_pts_2D(2,:),X_pts_2D(1,:),'+r');
		end
		% Calculate if the light source is on the same side
		%	of the plane to the optical center
		same_side = transpose(Xs(:,i_p)-S)*Xs(:,i_p);
		diff = diff + 10*(abs(same_side)-same_side);
		% Copy the final results error in the general vectors
		err((1+(i_p-1)*n_iso*nb_pts*2):(i_p*n_iso*nb_pts*2)) = diff(:);
	end

end


function [err] = vectorized_function_opt_d_parameters(K,n_iso,n_p,pts_2D,nb_pts,param_vect)
	S = param_vect(1:3);
	Rs = eul2rotm(transpose(reshape(param_vect(4:(4+(3*n_p-1))),3,n_p)));
	d_vector = transpose(param_vect((4+(3*n_p)):(3+4*n_p)));
	h_vector = -d_vector-dot(S*ones(1,n_p),reshape(Rs(:,3,:),3,n_p));
	r_vector = reshape(param_vect((4+4*n_p):(3+(n_iso+4)*n_p)),n_iso,n_p);
	Xs = S+h_vector.*reshape(Rs(:,3,:),3,n_p);
	angles = reshape(param_vect((4+(n_iso+4)*n_p):end),nb_pts*n_iso,n_p);
	cos_angles = cos(angles);
	sin_angles = sin(angles);
	zeros_angles = zeros(size(angles));
	err = zeros(n_iso*n_p*2*nb_pts,1);
	for i_p=1:n_p
		r_angles = repmat(transpose(r_vector(:,i_p)),nb_pts,1);
		X_pts = [r_angles(:).*cos_angles(:,i_p),r_angles(:).*sin_angles(:,i_p),zeros_angles(:,i_p)]*transpose(Rs(:,:,i_p));
		X_pts_2D = K*(Xs(:,i_p)+transpose(X_pts));
		X_pts_2D = X_pts_2D([2,1],:)./X_pts_2D(3,:);
		diff = pts_2D((1+(i_p-1)*n_iso*nb_pts):i_p*n_iso*nb_pts,:)-transpose(X_pts_2D);
		% Calculate if the light source is on the same side
		%	of the plane to the optical center
		same_side = transpose(Xs(:,i_p)-S)*Xs(:,i_p);
		diff = diff + 10*(abs(same_side)-same_side);
		% Copy the final results error in the general vectors
		err((1+(i_p-1)*n_iso*nb_pts*2):(i_p*n_iso*nb_pts*2)) = diff(:);
	end

end

%function [err] = vectorized_function_opt_source_fixed(S,K,n_iso,n_p,pts_2D,nb_pts,param_vect)
%	Rs = eul2rotm(transpose(reshape(param_vect(1:(3*n_p)),3,n_p)));
%	h_vector = transpose(param_vect((1+(3*n_p)):(4*n_p)));
%	r_vector = reshape(param_vect((1+4*n_p):((n_iso+4)*n_p)),n_iso,n_p);
%	Xs = S+h_vector.*reshape(Rs(:,3,:),3,n_p);
%	angles = reshape(param_vect((1+(n_iso+4)*n_p):end),nb_pts*n_iso,n_p);
%	cos_angles = cos(angles);
%	sin_angles = sin(angles);
%	zeros_angles = zeros(size(angles));
%	err = zeros(n_iso*n_p*2*nb_pts,1);
%	for i_p=1:n_p
%		r_angles = repmat(transpose(r_vector(:,i_p)),nb_pts,1);
%		X_pts = [r_angles(:).*cos_angles(:,i_p),r_angles(:).*sin_angles(:,i_p),zeros_angles(:,i_p)]*transpose(Rs(:,:,i_p));
%		X_pts_2D = K*(Xs(:,i_p)+transpose(X_pts));
%		X_pts_2D = X_pts_2D([2,1],:)./X_pts_2D(3,:);
%		diff = pts_2D((1+(i_p-1)*n_iso*nb_pts):i_p*n_iso*nb_pts,:)-transpose(X_pts_2D);
%		err((1+(i_p-1)*n_iso*nb_pts*2):(i_p*n_iso*nb_pts*2)) = diff(:);
%	end
%end


function [cell_arg] = param_vect_to_cell(param_vect,K,n_iso,n_p,data)
	cell_arg = cell(1,5+2*n_p*n_iso);
        cell_arg{1} = data.K;
        cell_arg{2} = param_vect(1:3);
	i_current = 4;
        cell_arg{3} = reshape(param_vect(i_current:(i_current+(3*n_p-1))),3,n_p);
	i_current = i_current+(3*n_p);
        cell_arg{4} = reshape(param_vect(i_current:(i_current+(n_p-1))),1,n_p);
	i_current = i_current+n_p;
        cell_arg{5} = reshape(param_vect(i_current:(i_current+n_iso*n_p-1)),n_iso,n_p);
	i_current = i_current+n_iso*n_p;
	i_cur = 6;
        for  i_p=1:n_p
                for i_iso = 1:n_iso
                        cell_arg{i_cur} = data.isocontour.Points{i_p}{i_iso};
			nb_pt = length(data.isocontour.Points{i_p}{i_iso});
                        cell_arg{i_cur+1} = reshape(param_vect(i_current:(i_current+nb_pt-1)),1,nb_pt);
			i_cur = i_cur+2;
			i_current = i_current+nb_pt;
                end
        end
end

function [param_vect,J_sparse] = cell_to_param_vect(varargin)
        K = varargin{1};
        S = varargin{2};
	param_vect = [];
	param_vect = [param_vect;S(:)];
        angles_euler = varargin{3};
	param_vect = [param_vect;angles_euler(:)];
        h_vector = varargin{4};
	param_vect = [param_vect;h_vector(:)];
        n_plane = size(angles_euler,2);
        r_vector = varargin{5};
	param_vect = [param_vect;r_vector(:)];
        m_iso = size(r_vector,1);
	nb_err = 0;
        for n=1:n_plane
		for m=1:m_iso
			i_current = 6+(n-1)*m_iso*2+(m-1)*2;
			angles = varargin{i_current+1};
			param_vect = [param_vect;angles(:)];
			nb_err = nb_err+2*length(angles);
                end
        end
	if nargout > 1
		nb_param = 3+(m_iso+4)*n_plane+(nb_err/2);
		i_0_param = 3+(m_iso+4)*n_plane;
		J_sparse = sparse(nb_err,nb_param);
		J_sparse(:,1:(3+(m_iso+4)*n_plane)) = 1;
		for i = 1:nb_err
			J_sparse(i,ceil(i/2)+i_0_param) = 1;
		end
	end
end
