% This function will be used to estimate the pose parameters directly with dense photometric approach
% Input Parameters :
%	* data : struct containing information like intrinsic parameters and groudntruth
%	* pts_cell : cell of set of points in the image for each surface
%	* I_pt_cell : cell of set of intensity values corresponding to the points in pts_cell
%	* dssp_init : initialization for the distance of the plane to the source
%	* psp_init : initialization for the pose of the planes
%	* ps_init : initialization for the position of the light source
function [dssp_opt,psp_opt,ps_opt,err_out] = global_dense_isocontours_optimization(...
	data,pts_cell,I_pt_cell,I_vect_cell,dssp_init,psp_init,ps_init)

	NB_R = length(I_vect_cell{1});
	nb_plane = length(pts_cell);
	% Parameters to optimize :
	% 	* Planes Normal, Ns
	%	* Planes distance to origin, ds
	%	* Light Source position, S
	%	* Spline control points (radius only, intensity can be fixed) rss
	% Initialization des variables
	% Vectoriser les données d'entrées
	inv_K = inv(data.K);
	% Getting initialize parameters
	%	Check if only one single intensity vector for the spline is given
	%	In that case use the single spline model
	if length(I_vect_cell)<length(pts_cell)
		single_spline = true;
	else
		single_spline = false;
	end
	[x_init,J_sparse,H_sparse,constraints_vector] = cell_to_init_param_vect(...
		inv_K,ps_init,psp_init,pts_cell,NB_R,single_spline);
	% Optimization setup for least square problem
	%options = optimoptions(@lsqnonlin,...
	%	'Algorithm','levenberg-marquardt',...
	%	'SpecifyObjectiveGradient',false,...
	%	'CheckGradients',false,...
	%	'FiniteDifferenceType','central',...
	%	'FiniteDifferenceStepSize',1e-7,...
	%	'UseParallel',true,...
	%	'Display','iter-detailed',...
	%	'JacobPattern',J_err_global);
	%fun_opt = @(x_param) evaluate_error_global_dense_derivate(inv_K,pts_cell,I_pt_cell,I_vect_cell,x_param,0);
	%tic;
	%x_opt = lsqnonlin(fun_opt,x_init,[],[],options);
	%toc
	%evaluate_error_global_dense_derivate(inv_K,pts_cell,I_pt_cell,I_vect_cell,x_init,1,single_spline);

	% Optimization setup for norm 2 or norm 1 direct problem
	options = optimoptions(@fmincon,...
		'SpecifyObjectiveGradient',true,...
		'CheckGradients',false,...
		'FiniteDifferenceStepSize',1e-6,...
		'MaxIterations',400,...
		'FiniteDifferenceType','central',...
		'UseParallel',true,...
		'Display','off',...
		'HessPattern',H_sparse);
	%'FiniteDifferenceType','central',...
	fun_opt = @(x_param) evaluate_norm_1_global_dense_derivate(inv_K,pts_cell,...
		I_pt_cell,I_vect_cell,x_param,single_spline);
	%evaluate_error_global_dense_derivate(inv_K,pts_cell,I_pt_cell,I_vect_cell,x_init,1);


	[err_global] = evaluate_error_global_dense_derivate(inv_K,pts_cell,I_pt_cell,I_vect_cell,x_init,1,single_spline);
	tic;
	% fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options)
	x_opt = fmincon(fun_opt,x_init,[],[],[],[],constraints_vector,[],[],options);
	toc;
	%x_opt = x_init;
	%dssp = dssp_init;
	%psp = psp_init;
	%ps = ps_init;
	[err_global] = evaluate_error_global_dense_derivate(inv_K,pts_cell,I_pt_cell,I_vect_cell,x_opt,0,single_spline);
	err_out = sum(err_global.^2);
	%fun_opt = @(x) norm(evaluate_error_global_dense_derivate(inv_K,pts_cell,I_pt_cell,I_vect_cell,x,0));
	%options = optimoptions('fminunc','Algorithm',...
	%	'quasi-newton','Display','iter-detailed');
	%%x_opt = fminunc(fun_opt,x_init,options);
	%x_opt = x_init;
	%tic;
	[S,psp_opt,r_vects] = param_vect_to_cell(x_opt,nb_plane,NB_R,single_spline);
	ps_opt = {S};
	dssp_opt = calculate_dssp(ps_opt,psp_opt);
end

function [dssp] = calculate_dssp(ps,psp)
	nb_plane = length(psp);
	dssp = cell(1,nb_plane);
	for i_p = 1:nb_plane
		dssp{i_p} = {-(psp{i_p}{1}(4)+transpose(psp{i_p}{1}(1:3))*ps{1})};
	end
end

function [S,psp,r_vects] = param_vect_to_cell(param_vect,nb_plane,NB_R,single_spline)
	S = param_vect(1:3);
	psp = cell(1,nb_plane);
	r_vects = cell(1,nb_plane);
	J_r_vect = tril(ones(NB_R));
	if nargin<4 || ~single_spline
		single_spline = false;
		i_param = 4;
	else
		i_param = 4+NB_R;
		a_vect = transpose(J_r_vect)*param_vect(4:(3+NB_R));
	end
	for i_p	= 1:nb_plane
		if ~single_spline
			N = angle_to_normal(param_vect(i_param),param_vect(i_param+1));
			d = param_vect(i_param+2);
			psp{i_p} = {[N;d]};
			r_vects{i_p} = transpose(J_r_vect*param_vect((i_param+3):(i_param+2+NB_R)));
			i_param = i_param+3+NB_R;
		else
			N = angle_to_normal(param_vect(i_param),param_vect(i_param+1));
			d = param_vect(i_param+2);
			psp{i_p} = {[N;d]};
			r_vects{i_p} = a_vect;
			i_param = i_param+3;
		end
	end
end

function [norm_L1,gradient_norm_L1] = evaluate_norm_1_global_dense_derivate(inv_K,pts_cell,I_pt_cell,I_vect_cell,x_param,single_spline)
	if nargin < 6
		single_spline = false;
	end
	if nargout > 1
		[err,J_err] = evaluate_error_global_dense_derivate(inv_K,pts_cell,I_pt_cell,I_vect_cell,x_param,0,single_spline);
		[norm_L1,gradient_norm_L1] = error_norm_L1(err,J_err);
	else
		err = evaluate_error_global_dense_derivate(inv_K,pts_cell,I_pt_cell,I_vect_cell,x_param,0,single_spline);
		[norm_L1] = error_norm_L1(err);
	end
end
