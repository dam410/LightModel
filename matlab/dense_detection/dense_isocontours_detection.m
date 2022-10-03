% Function to estimate the homography for each plane individually based on our isocontour model
%	using splien interpolation, the optimization is dense on each pixels of the image representing the surface
% Input Parameters :
%	I : Image undistorted and filtered
%	K : Intrinsic camera parameter matrix
%	polys_2D : Cell of polygons delimitation of each surface in the image
% Optional input parameters
%	mode : String parametering the type of optimization process we apply
%	display : String ("On"/"Off")  Debug and displaying the curves
%	homographies : Cell of the initial homography
% Output Parameters :
%	homogs_1 : Cell of the output homographies
%	homogs_2 : Cell of the second output homographies if we use ambiguity (Same as homogs_1 otherwise)
%	spline_params : Cell of the control points of the spline for each plane
function [homogs_1,homogs_2,spline_params,optim_results] = dense_isocontours_detection(I,K,polys_2D,options)
	arguments
                I double
                K (3,3) double
		polys_2D (:,:) cell
		options.Mode (1,1) string = 'perspective'
		options.Display (1,1) string = 'off'
		options.Nb_R (1,1) double = 5;
		options.Homographies (:,:) cell
		options.Points (:,:) cell
		options.Intensities (:,:) cell
		options.N_circ (1,1) double = 3;
        end
	if size(I,3)==3
                I = double(rgb2gray(uint8(I)));
        else
                I = double(I);
        end
	[n,m] = size(I);
	n_p = length(polys_2D);
	%I_gaussian = I;
	I = wiener2(I,[15,15]);
	I_gaussian = imgaussfilt(I,3,'FilterSize',11);
	% Mode to parameter the type of optimization we want
	mode = options.Mode;
	display = options.Display;
	% Number of radius to interpolate the spline
	NB_R = options.Nb_R;
	% Number of extracted isocontours
	N_CIRC = options.N_circ;
	% Inverse matrix to pass from [r0,r1,r2] -> [r0,a1,a2]
	%	r1 = r0+a1 ... r2 = r0+a1+a2
	% This allows to add simple constraints [r0,...,an] >= [0,..,0]
	R_vect = tril(ones(NB_R));
	inv_R_vect = inv(R_vect);
	% Return parameters the homographies with ambiguity if asked
	ambig = 'no';
	homogs_1 = cell(1,n_p);
	if nargout > 1
		ambig = 'yes';
		homogs_2 = cell(1,n_p);
		if nargout > 2
			spline_params = cell(1,n_p);
		end
	end
	optim_results = cell(1,n_p);
	pt_barycenters = cell(1,n_p);
	pt_centerss = cell(1,n_p);
	is_data_normalized = true;
	switch mode
		case 'perspective'
			initialization_type = 'rough';
			first_opti = 'perspective';
			ippe = 'no';
			second_opti = 'no';
		case 'perspective_L1'
			initialization_type = 'groundtruth_inverse';
			first_opti = 'perspective_L1';
			ippe = 'no';
			second_opti = 'no';
		case 'perspective_L2'
			initialization_type = 'rough';
			first_opti = 'perspective_L2';
			ippe = 'no';
			second_opti = 'no';
		case 'perspectivity'
			initialization_type = 'rough';
			first_opti = 'perspectivity';
			ippe = 'no';
			%second_opti = 'perspectivity';
			second_opti = 'no';
		case 'rotation'
			initialization_type = 'rough';
			first_opti = 'rotation';
			ippe = 'no';
			is_data_normalized = false;
			%second_opti = 'perspectivity';
			second_opti = 'no';
		case 'perspectivity_no_ambiguity'
			initialization_type = 'rough';
			first_opti = 'perspectivity';
			ippe = 'no';
			second_opti = 'no_second_opt';
			%second_opti = 'no';
		case 'perspectivity_L1'
			initialization_type = 'rough';
			first_opti = 'perspectivity_L1';
			ippe = 'no';
			second_opti = 'no';
		case 'affine'
			initialization_type = 'rough';
			first_opti = 'affine';
			ippe = 'no';
			second_opti = 'no';
		case 'ippe_no_opt'
			%initialization_type = 'groundtruth';
			initialization_type = 'rough';
			%first_opti = 'affine';
			first_opti = 'perspective';
			ippe = 'no';
			second_opti = 'no';
		case 'ippe_perspective'
			initialization_type = 'rough';
			first_opti = 'affine';
			ippe = 'yes';
			second_opti = 'perspective';
		case 'groundtruth_perspective'
			initialization_type = 'groundtruth_inverse';
			first_opti = 'perspective';
			ippe = 'no';
			second_opti = 'no';
		case 'groundtruth_translation'
			initialization_type = 'groundtruth';
			first_opti = 'translation'
			ippe = 'no';
			second_opti = 'no';
		case 'groundtruth_inverse_translation'
			initialization_type = 'groundtruth_inverse';
			first_opti = 'translation'
			ippe = 'no';
			second_opti = 'no';
		otherwise
			initialization_type = 'rough';
			first_opti = 'perspective';
			ippe = 'no';
			second_opti = 'no';
	end

	% Calculate for each plane the best isocontour model
	for i_p = 1:length(polys_2D)
		if ~(isfield(options,'Points') && isfield(options,'Intensities'))
			[pt_in_poly,I_pt,pt_barycenter,I_vect,ind_pt] = img_points_from_poly(I_gaussian,polys_2D{i_p},NB_R);
			pt_barycenters{i_p} = pt_barycenter;
			pt_barycenter = [pt_barycenter(1),pt_barycenter(2)];
		else
			disp('Using provided intensities and points');
			pt_in_poly = options.Points{i_p};
			I_pt = options.Intensities{i_p};
			% Select the intensity range
			[I_max,ind_max] = max(I_pt);
			[I_min,ind_min] = min(I_pt);
			I_vect = I_max:-(I_max-I_min)/(NB_R-1):I_min;
			% Barycenter and point index (For display purpose)
			pt_barycenter = mean(pt_in_poly);
			ind_pt = sub2ind([n,m],pt_in_poly(:,2),pt_in_poly(:,1));
		end
		if strcmp(display,'on');
			I_test = I_gaussian;
			I_test(ind_pt) = 0;
			figure('Name','Masked area');
			imshow(I_test/255);
		end

		%% Data normalization if needed, not used in perspectivity or rotation case H_normalized = eye(3)
		if is_data_normalized
			pt_in_centered = pt_in_poly-pt_barycenter;
			var_pt = var(pt_in_centered);
			pt_in_normalized = [pt_in_centered(:,1)./var_pt(1),pt_in_centered(:,2)./var_pt(2)];
			H_normalized = [var_pt(1),0,pt_barycenter(1);0,var_pt(2),pt_barycenter(2);0,0,1];
		else
			pt_in_normalized = pt_in_poly;
			H_normalized = eye(3);

		end

		% Calculate the initial parameters
		% 	If nargin>3, homography is given, we are using it as initialization
		%		perpective homography
		%	Otherwise, we use initialization with pt_min and pt_max
		%		affine homography
		switch initialization_type
		case 'groundtruth'
			if ~isfield(options,'Homographies')
				error('Homographies not given for groundtruth initialization');
			end
			[h_vect_init,r_vect_init] = initial_parameters_given_homography(...
				I,pt_in_normalized,I_vect,inv(H_normalized)*options.Homographies{i_p});
		case 'groundtruth_inverse'
			if ~isfield(options,'Homographies')
				error('Homographies not given for groundtruth initialization');
			end
			[h_vect_init,r_vect_init] = initial_parameters_given_homography_inverse(...
				I,pt_in_normalized,I_vect,options.Homographies{i_p}*H_normalized);
		case 'rough'
			[h_vect_init,r_vect_init] = initial_parameters_rough(...
				I,pt_in_poly,pt_in_normalized,I_vect,I_pt);
		end
		a_vect_init = inv_R_vect*r_vect_init;

		% Calculate and display the variance of the points before the optimisation
		[I_pt_proj,r_proj_init,err_opt] = evaluate_error_homography_monotone_detailed(...
                        I,pt_in_normalized,I_pt,I_vect,vector_to_homography_perspective(h_vect_init),...
			inv_R_vect*r_vect_init);
                %display_spline_function(r_vect_init,r_proj_init,I_vect,I_pt,...
		%	'Display the variance before optimisation');


		%% Display initial error in 3D
		if strcmp(display,'on');
			[I_pt_init,r_proj_init,err_init] = ...
				evaluate_error_homography_monotone_detailed(I,pt_in_normalized,I_pt,I_vect,...
				vector_to_homography_perspective(h_vect_init),a_vect_init);
			display_error_3D(err_init,ind_pt,polys_2D{i_p},I,...
					'Difference between image and our interpolation initially');
			display_spline_function(r_vect_init,r_proj_init,I_vect,I_pt,...
					'Display spline interpolation initialization');	
			pause(0.5);
		end


		% Options parameters for the optimization
		options_opt = optimoptions('lsqnonlin');
		options_opt.Algorithm = 'levenberg-marquardt';
		options_opt.Display = 'iter-detailed';
		options_opt.StepTolerance = 1e-12;
		options_opt.FunctionTolerance = 1e-10;
		options_opt.MaxIter = 80;

		% First optimization affine or perspective
		% Calculate the optimal homography H_opt and optimal radius r_opt
		switch first_opti 
		case 'perspective'
			options_perspective = optimoptions('lsqnonlin');
			options_perspective.Algorithm = 'trust-region-reflective';
			options_perspective.Display = 'iter-detailed';
			options_perspective.StepTolerance = 1e-12;
			options_perspective.FunctionTolerance = 1e-10;
			options_perspective.MaxIter = 200;
			options_perspective.SpecifyObjectiveGradient = false;
			options_perspective.CheckGradients = false;
			options_perspective.FiniteDifferenceType = 'central';
			options_perspective.FiniteDifferenceStepSize = 1e-8;
			fun_opt = @(x_param) fun_opt_persp(I,pt_in_normalized,I_pt,I_vect,x_param);
			x_init = [h_vect_init;a_vect_init];
			disp(['First optimization with ',first_opti,' homography']);
			% Define homography upper and lower bound
			lb = [-Inf*ones(8,1);[-Inf;zeros(NB_R-1,1)]];
			ub = Inf*ones(8+NB_R,1);
			tic;
			try
				[x_opt,~,~,~,output] = lsqnonlin(fun_opt,x_init,lb,ub,options_perspective);
				%[x_opt,~,~,~,output] = lsqnonlin(fun_opt,x_init,[],[],options_perspective);
			catch exception
				x_opt = x_init;
				output = struct();
				output.message = 'Failed to optimized LSQ Pb';
			end
			toc
			H_opt_normalized = vector_to_homography_perspective(x_opt(1:8));
			r_opt = R_vect*x_opt(9:end);
			optim_results{i_p} = output;
		case 'perspective_L1'
			%%%%% First optimization only with geometry
			%options_L1 = optimoptions('fminunc');
			%options_L1.SpecifyObjectiveGradient = true;
			%options_L1.Display = 'iter-detailed';
			%options_L1.FunctionTolerance = 1e-10;
			%options_L1.StepTolerance = 1e-12;
			%options_L1.FiniteDifferenceStepSize = 1e-11;
			%options_L1.MaxIter = 20;
			%options_L1.CheckGradients = false;
			%options_L1.Algorithm = 'quasi-newton';
			%fun_opt = @(param) fun_opt_persp_L1_only_h(I,pt_in_normalized,I_pt,I_vect,[param;a_vect_init]);
			%disp(['First optimization with L1 norm and perspective homography']);
			%tic;
			%x_opt = fminunc(fun_opt,h_vect_init,options_L1);
			%toc
			%%%% Second optimization with photometry and geometry
			options_L1 = optimoptions('fmincon');
			options_L1.SpecifyObjectiveGradient = true;
			options_L1.Display = 'iter-detailed';
			options_L1.FunctionTolerance = 1e-10;
			options_L1.StepTolerance = 1e-12;
			options_L1.FiniteDifferenceStepSize = 1e-11;
			options_L1.FiniteDifferenceType = 'central';
                        options_L1.FiniteDifferenceStepSize = 1e-8;
			options_L1.MaxIter = 200;
			options_L1.CheckGradients = false;
			options_L1.Algorithm = 'interior-point';
			lb = [-Inf*ones(8,1);[-Inf;zeros(NB_R-1,1)]];
			ub = Inf*ones(8+NB_R,1);
			fun_opt = @(param) fun_opt_persp_L1(I,pt_in_normalized,I_pt,I_vect,param);
			%x_init = [x_opt(1:8);a_vect_init];
			x_init = [h_vect_init;a_vect_init];
			tic;
			[x_opt,~,~,output] = fmincon(fun_opt,x_init,[],[],[],[],lb,ub,[],options_L1);
			toc
			H_opt_normalized = vector_to_homography_perspective(x_opt(1:8));
			r_opt = R_vect*x_opt(9:end);
			optim_results{i_p} = output;
		case 'perspective_L2'
			%%%%% First optimization only with geometry
			%options_L2 = optimoptions('fminunc');
			%options_L2.SpecifyObjectiveGradient = true;
			%options_L2.Display = 'iter-detailed';
			%options_L2.FunctionTolerance = 1e-10;
			%options_L2.StepTolerance = 1e-12;
			%options_L2.FiniteDifferenceStepSize = 1e-11;
			%options_L2.MaxIter = 20;
			%options_L2.CheckGradients = false;
			%options_L2.Algorithm = 'quasi-newton';
			%fun_opt = @(param) fun_opt_persp_L2_only_h(I,pt_in_normalized,I_pt,I_vect,[param;a_vect_init]);
			%disp(['First optimization with L2 norm and perspective homography']);
			%tic;
			%x_opt = fminunc(fun_opt,h_vect_init,options_L2);
			%toc
			%%%% Second optimization with photometry and geometry
			options_L2 = optimoptions('fmincon');
			options_L2.SpecifyObjectiveGradient = true;
			options_L2.Display = 'iter-detailed';
			options_L2.FunctionTolerance = 1e-10;
			options_L2.StepTolerance = 1e-12;
			options_L2.FiniteDifferenceStepSize = 1e-11;
			options_L2.MaxIter = 200;
			options_L2.FiniteDifferenceType = 'central';
			options_L2.FiniteDifferenceStepSize = 1e-8;
			options_L2.CheckGradients = false;
			options_L2.Algorithm = 'interior-point';
			lb = [-Inf*ones(8,1);[-Inf;zeros(NB_R-1,1)]];
			ub = Inf*ones(8+NB_R,1);
			fun_opt = @(param) fun_opt_persp_L2(I,pt_in_normalized,I_pt,I_vect,param);
			%x_init = [x_opt(1:8);a_vect_init];
			x_init = [h_vect_init;a_vect_init];
			tic;
			[x_opt,~,~,output] = fmincon(fun_opt,x_init,[],[],[],[],lb,ub,[],options_L2);
			toc
			H_opt_normalized = vector_to_homography_perspective(x_opt(1:8));
			r_opt = R_vect*x_opt(9:end);
			optim_results{i_p} = output;
		case 'perspectivity'
			h_vect_init_not_normalized = homography_to_vector_perspective(...
				vector_to_homography_perspective(h_vect_init)*inv(H_normalized));
			[v_init_persp,scale] = vector_init_to_perspectivity_minimal_param(h_vect_init_not_normalized,inv(K));
			x_init = [v_init_persp;a_vect_init./(scale.^2)];
			% Method without changing intensity in control point
			[H_opt_normalized,r_opt] = optimize_perspectivity_monotonic(I,pt_in_poly,I_pt,I_vect,...
				inv(K),x_init,H_normalized,NB_R,R_vect);
			% Changing intensity in control point
			%[H_opt_normalized,r_opt,I_vect_opt] = optimize_perspectivity_full_cp(I,pt_in_poly,I_pt,I_vect,...
			%	inv(K),x_init,H_normalized,NB_R,R_vect);
		case 'rotation'
			% Model based on simple rotation of the main axis H = R
			% Note : There is no normalization in this case H_normalized = eye(3)
			h_vect_init_not_normalized = homography_to_vector_perspective(...
                                vector_to_homography_perspective(h_vect_init)*inv(H_normalized));
			[v_init_persp,scale] = vector_init_to_perspectivity_minimal_param(...
				h_vect_init_not_normalized,inv(K));
                        x_init = [v_init_persp;a_vect_init./(scale.^2)];
                        [v_init_persp,scale] = vector_init_to_rotation_param(h_vect_init,inv(K));
                        x_init = [v_init_persp;a_vect_init./(scale.^2)];
                        [H_opt_normalized,r_opt] = optimize_rotation_only(I,pt_in_poly,I_pt,I_vect,...
                                inv(K),x_init,NB_R,R_vect);
		case 'perspectivity_L1'
			options_L1 = optimoptions('fmincon');
			options_L1.SpecifyObjectiveGradient = false;
			options_L1.Display = 'iter-detailed';
			options_L1.FunctionTolerance = 1e-10;
			options_L1.StepTolerance = 1e-12;
			options_L1.FiniteDifferenceStepSize = 1e-11;
			options_L1.MaxIter = 200;
			options_L1.FiniteDifferenceType = 'central';
			options_L1.FiniteDifferenceStepSize = 1e-8;
			options_L1.CheckGradients = false;
			options_L1.Algorithm = 'interior-point';
			lb = [-Inf;0;ones(2,1);[-Inf;zeros(NB_R-1,1)]];
			ub = [Inf;pi;ones(2+NB_R,1)];
			invK = inv(K);
                        fun_opt = @(param) error_norm_L1(evaluate_error_homography_monotone(I,pt_in_poly,I_pt,I_vect,...
                                vector_to_perspectivity_minimal_param(param(1:4),invK),param(5:(4+NB_R))));
                        h_vect_init_not_normalized = homography_to_vector_perspective(...
                                vector_to_homography_perspective(h_vect_init)*inv(H_normalized));
                        [v_init_persp,scale] = vector_init_to_perspectivity_minimal_param(h_vect_init_not_normalized,invK);
                        x_init = [v_init_persp;a_vect_init./scale.^2];
			tic;
			[x_opt,~,~,output] = fmincon(fun_opt,x_init,[],[],[],[],lb,ub,[],options_L1);
			toc
			H_opt_normalized = vector_to_perspectivity_minimal_param(x_opt(1:4),invK)*H_normalized;
                        r_opt = R_vect*x_opt(5:(4+NB_R));
			optim_results{i_p} = output;
		case 'affine'
			options_affine = options_opt;
			options_affine.SpecifyObjectiveGradient = true;
			fun_opt = @(x_param) fun_opt_affine(I,pt_in_normalized,I_pt,I_vect,x_param);
			x_init = [homography_to_vector_affine(...
				vector_to_homography_perspective(h_vect_init));a_vect_init];
			homography_to_vector_affine(...
				vector_to_homography_perspective(h_vect_init))
			disp(['First optimization with ',first_opti,' homography']);
			tic;
			x_opt = lsqnonlin(fun_opt,x_init,[],[],options_affine);
			toc
			H_opt_normalized = vector_to_homography_affine(x_opt(1:6));
			r_opt = R_vect*x_opt(7:end);
		case 'translation'
			options_opt = optimoptions('lsqnonlin');
			options_opt.Algorithm = 'levenberg-marquardt';
			options_opt.SpecifyObjectiveGradient = true;
			options_opt.Display = 'iter-detailed';
			options_opt.StepTolerance = 1e-12;
			options_opt.FunctionTolerance = 1e-10;
			options_opt.MaxIter = 80;
			H_init = vector_to_homography_perspective(h_vect_init);
			% TODO Calculate the new pt_in_poly so it will be easier to calculate 
			% derivative depending on p1 and p2
			fun_opt = @(param) fun_opt_translation(I,pt_in_normalized,I_pt,I_vect,H_init,param);
			%fun_opt = @(param) evaluate_error_homography_monotone(I,pt_in_normalized,I_pt,I_vect,...
			%	[1,0,param(1);0,1,param(2);0,0,1]*H_init,param(3:end));
			x_init = [0;0;a_vect_init];
			disp(['First optimization with ',first_opti,' homography']);
			tic;
			x_opt = lsqnonlin(fun_opt,x_init,[],[],options_opt);
			toc
			H_opt_normalized = [1,0,x_opt(1);0,1,x_opt(2);0,0,1]*H_init;
			r_opt = R_vect*x_opt(3:end);
		case 'photometric'
			H_init = vector_to_homography_perspective(h_vect_init);
			fun_opt = @(param) evaluate_error_homography_monotone(I,pt_in_normalized,I_pt,I_vect,...
				H_init,param(1:end));
			x_init = [a_vect_init];
			disp(['First optimization with ',first_opti,' homography']);
			tic;
			x_opt = lsqnonlin(fun_opt,x_init,[],[],options_opt);
			toc
			H_opt_normalized = H_init;
			r_opt = R_vect*x_opt;
		case 'no'
			x_opt = [h_vect_init;a_vect_init];
			H_opt_normalized = vector_to_homography_perspective(x_opt(1:8));
			r_opt = R_vect*x_opt(9:end);
		end
		H_opt = H_opt_normalized*inv(H_normalized);
		% Just reaffect all control points intensity with the optimized one
		if exist('I_vect_opt')
			I_vect = I_vect_opt;
		end

		% Calculate the error after first optimization
		[I_pt_proj,r_proj,err_opt] = evaluate_error_homography_monotone_detailed(...
			I,pt_in_poly,I_pt,I_vect,H_opt,inv_R_vect*r_opt);
		%display_spline_function(r_opt,r_proj,I_vect,I_pt,'Display the variance after optimisation');


		% Displaying error
		if ~strcmp(first_opti,'no') && strcmp(display,'on')
			% Display error in 3D after first optimization
			display_error_3D(err_opt,ind_pt,polys_2D{i_p},I,...
			'Difference between image and our interpolation after first optimization');
			pause(0.5);

			% Display the spline interpolated
			display_spline_function(r_opt,r_proj,I_vect,I_pt,...
			'Spline function interpolated after first optimization');
			
			pause(0.5);
		end

		% IPPE algorithme is used to calculate the two ambiguous 3D perspective plane
		%	poses that can be obtained from the affine first optimization	
		% This is justified by the fact that if the surface is small and far from the camera
		%	the homography can be ambiguous, this can make the convergence basin revolves
		%	around two local minimum
		% To solve that we calculate two initial point that are perspectively opposite and 
		%	start the optimization from both positions
		switch ippe
		case 'yes'
			%[h_vects,pt_centers] = homography_to_perspectivity(K,H_opt,pt_barycenters{i_p});
			invK = inv(K);
			%save('temp_debug_ippe.mat','invK','H_opt');
			h_vects = homography_to_ambiguous_perspectivity(invK,H_opt);
			switch second_opti
			case 'perspectivity'
			

				% 
				a_vect_init = inv_R_vect*r_opt;
				%[H_opt_normalized,r_opt] = optimize_perspectivity(I,pt_in_poly,I_pt,I_vect,...
				%	inv(K),x_init,H_normalized,NB_R,R_vect);
	
				% First ambiguity
				[v_init_persp_1,scale_1] = vector_init_to_perspectivity_minimal_param(h_vects(:,1),invK);
				x_init_1 = [v_init_persp_1;a_vect_init./(scale_1.^2)];
				[H_opt_1,r_opt_1] = optimize_perspectivity(I,pt_in_poly,I_pt,I_vect,...
					invK,x_init_1,eye(3),NB_R,R_vect);
				% Second ambiguity
				[v_init_persp_2,scale_2] = vector_init_to_perspectivity_minimal_param(h_vects(:,2),invK);
				x_init_2 = [v_init_persp_2;a_vect_init./(scale_2.^2)];
				[H_opt_2,r_opt_2] = optimize_perspectivity(I,pt_in_poly,I_pt,I_vect,...
					invK,x_init_2,eye(3),NB_R,R_vect);
				% Put them all together
				h_vects_opt = [homography_to_vector_perspective(H_opt_1),...
					homography_to_vector_perspective(H_opt_2)];
				r_vects_opt = [r_opt_1,r_opt_2];

			case 'no'
				% No optimization
				h_vects_opt = h_vects;
				r_vects_opt = [r_opt,r_opt];
			end
			switch ambig
			case 'no'
				% Keep the best position
				err_12 = zeros(1,2);
				err_12(1) = norm(evaluate_error_homography_monotone(I,pt_in_poly,I_pt,...
					I_vect,vector_to_homography_perspective(h_vects_opt(:,1)),...
					inv_R_vect*r_vects_opt(:,1)));
				err_12(2) = norm(evaluate_error_homography_monotone(I,pt_in_poly,I_pt,...
					I_vect,vector_to_homography_perspective(h_vects_opt(:,2)),...
					inv_R_vect*r_vects_opt(:,2)));
				[err_norm,best_12] = min(err_12);
				% In case of failure to use the IPPE algorithm
				if (err_norm > err_opt) & strcmp('perspective',second_opti)
					x_init = [homography_to_vector_perspective(H_opt);inv_R_vect*r_opt];
					tic;
					x_opt = lsqnonlin(fun_opt,x_init,[],[],options_opt);
					toc
					H_opt = vector_to_homography_perspective(x_opt(1:8));
					r_opt = R_vect*x_opt(9:end);				
				else
					r_opt = r_vects_opt(:,best_12);
					H_opt = vector_to_homography_perspective(h_vects(:,best_12));
				end
			case 'yes'
				r_opt = r_vects_opt(:,1);
				H_opt = vector_to_homography_perspective(h_vects(:,1));
				r_opt_2 = r_vects_opt(:,2);
				H_opt_2 = vector_to_homography_perspective(h_vects(:,2));
			end
		case 'no'
			if strcmp(ambig,'yes')
				r_opt_2 = r_opt;
				H_opt_2 = H_opt;
			end
		end

		% Calculate error and interpolated data
		[I_pt_proj,r_proj,err_opt] = evaluate_error_homography_monotone_detailed(...
			I,pt_in_poly,I_pt,I_vect,H_opt,inv_R_vect*r_opt);
		homogs_1{i_p} = H_opt;
		if strcmp(ambig,'yes')
			[I_pt_proj_2,r_proj_2,err_opt_2] = evaluate_error_homography_monotone_detailed(...
				I,pt_in_poly,I_pt,I_vect,H_opt_2,inv_R_vect*r_opt_2);
			homogs_2{i_p} = H_opt_2;
		end


		% Displaying error
		if (~strcmp(ippe,'no') || ~strcmp(second_opti,'no')) && strcmp(display,'on')
			% Display the spline interpolated
			display_spline_function(r_opt,r_proj,I_vect,I_pt,...
			'Spline function interpolated after final optimization');
			pause(0.5);

			% Display debug for error function
			display_error_3D(err_opt,ind_pt,polys_2D{i_p},I,...
			'Difference between image and our interpolation after final optimization');
			pause(0.5);

			if strcmp(ambig,'yes')
				% Display the spline interpolated
				display_spline_function(r_opt_2,r_proj_2,I_vect,I_pt,...
				'Spline function interpolated after final optimization');
				pause(0.5);

				% Display debug for error function
				display_error_3D(err_opt_2,ind_pt,polys_2D{i_p},I,...
				'Difference between image and our interpolation after final optimization');
				pause(0.5);	
			end
		end
		if nargout > 2
			spline_params{i_p} = [transpose(r_opt);I_vect];
		end	
	end
	
	% Display ellipses of selected radius
	if strcmp(display,'on')
		for i_p = 1:n_p
			% Calculate hypothetic points
			[points_iso,curves_iso] = isocontours_from_homography(...
				polys_2D{i_p},homogs_1{i_p},N_CIRC);
			points{i_p} = points_iso;
			curves{i_p} = curves_iso;
			if strcmp(ambig,'yes')
				[points_iso_2,curves_iso_2] = isocontours_from_homography(...
					polys_2D{i_p},homogs_2{i_p},N_CIRC);
				points_2{i_p} = points_iso_2;
				curves_2{i_p} = curves_iso_2;
			end
		end
		pause(0.8);
		figure('Name','Circles following homography');
		imshow(I/255);
		hold on;
		for i_p = 1:n_p
			for i_iso = 1:N_CIRC
				%plot(pt_barycenters{i_p}(1),pt_barycenters{i_p}(2),'+b');
				plot(points{i_p}{i_iso}(:,2),points{i_p}{i_iso}(:,1),'+r');
				if strcmp(ambig,'yes')
					plot(points_2{i_p}{i_iso}(:,2),points_2{i_p}{i_iso}(:,1),'-g');
					%if strcmp(ippe,'yes')
					%	plot(pt_centerss{i_p}{2}(1),pt_centerss{i_p}{2}(2),'+g');
					%end
				end
			end
		end
		pause(0.8);
	end
end

function display_error_3D(err,ind_pt,poly_2D,I,fig_name);
	[n,m] = size(I);
	min_x = max(1,min(m,floor(min(poly_2D(1,:)))));
        max_x = max(1,min(m,ceil(max(poly_2D(1,:)))));
        min_y = max(1,min(n,floor(min(poly_2D(2,:)))));
        max_y = max(1,min(n,ceil(max(poly_2D(2,:)))));

	err_matrix = zeros(size(I));
	err_matrix(ind_pt) = err;
	[X,Y] = meshgrid(1:m,1:n);
	X_short = X(min_y:max_y,min_x:max_x);	
	Y_short = Y(min_y:max_y,min_x:max_x);	
	err_matrix_short = err_matrix(min_y:max_y,min_x:max_x);	
	figure('Name',fig_name);
	surf(X_short,Y_short,err_matrix_short);
end

function [err,J_err] = fun_opt_affine(I,pt_in_poly,I_pt,I_vect,x_param)
        H_opt = vector_to_homography_affine(x_param(1:6));
        a_opt = x_param(7:end);
        nb_pt_interp = length(a_opt);
	if nargin>1
		[err,J_err_full] = evaluate_error_homography_monotone(...
			I,pt_in_poly,I_pt,I_vect,H_opt,a_opt);
		J_err = J_err_full(:,[1,2,4,5,7,8,7+(1:nb_pt_interp)]);
	else
		err = evaluate_error_homography_monotone(...
                        I,pt_in_poly,I_pt,I_vect,H_opt,a_opt)
	end
end

function [err,J_err] = fun_opt_translation(I,pt_in_poly,I_pt,I_vect,H_init,x_param)
	H_opt = [1,0,x_param(1);0,1,x_param(2);0,0,1]*H_init;
	a_opt = x_param(3:end);
	nb_pt_interp = length(a_opt);
	if nargin>1
		[err,J_err_full] = evaluate_error_homography_monotone(...
			I,pt_in_poly,I_pt,I_vect,H_opt,a_opt);
		% Special Jacobian matrix to only depends on the two translation parameters
		J_h_trans = [H_init(3,1),0;0,H_init(3,1);0,0;H_init(3,2),0;0,H_init(3,2);0,0;H_init(3,3),0;0,H_init(3,3);0,0];
		J_x_param = [J_h_trans,zeros(9,nb_pt_interp);zeros(nb_pt_interp,2),eye(nb_pt_interp)];
		J_err = J_err_full*J_x_param;
	else
		err = evaluate_error_homography_monotone(...
			I,pt_in_poly,I_pt,I_vect,H_opt,a_opt)
	end
end

function [err,J_err] = fun_opt_persp(I,pt_in_poly,I_pt,I_vect,x_param)
        H_opt = vector_to_homography_perspective(x_param(1:8));
        a_opt = x_param(9:end);
        nb_pt_interp = length(a_opt);
	if nargout>1
		[err,J_err_full] = evaluate_error_homography_monotone(...
			I,pt_in_poly,I_pt,I_vect,H_opt,a_opt);
		J_err = J_err_full(:,[1:8,9+(1:nb_pt_interp)]);
	else
		err = evaluate_error_homography_monotone(...
                        I,pt_in_poly,I_pt,I_vect,H_opt,a_opt);
	end
end


function [norm_L1,gradient_norm_L1] = fun_opt_persp_L1(I,pt_in_poly,I_pt,I_vect,x_param)
	if nargout > 1
		[err,J_err] = fun_opt_persp(I,pt_in_poly,I_pt,I_vect,x_param);
		[norm_L1,gradient_norm_L1] = error_norm_L1(err,J_err);
	else
		[err] = fun_opt_persp(I,pt_in_poly,I_pt,I_vect,x_param);
		[norm_L1] = error_norm_L1(err);
	end
end

function [norm_L1,gradient_norm_L1] = fun_opt_persp_L1_only_h(I,pt_in_poly,I_pt,I_vect,x_param)
	if nargout > 1
		[err,J_err] = fun_opt_persp(I,pt_in_poly,I_pt,I_vect,x_param);
		[norm_L1,gradient_norm_L1] = error_norm_L1(err,J_err);
		gradient_norm_L1 = gradient_norm_L1(1,1:8);
		%disp('Calling with J_err');
	else
		[err] = fun_opt_persp(I,pt_in_poly,I_pt,I_vect,x_param);
		[norm_L1] = error_norm_L1(err);
	end
end

function [norm_L2,gradient_norm_L2] = fun_opt_persp_L2(I,pt_in_poly,I_pt,I_vect,x_param)
	if nargout > 1
		[err,J_err] = fun_opt_persp(I,pt_in_poly,I_pt,I_vect,x_param);
		[norm_L2,gradient_norm_L2] = error_norm_L2(err,J_err);
	else
		[err] = fun_opt_persp(I,pt_in_poly,I_pt,I_vect,x_param);
		[norm_L2] = error_norm_L2(err);
	end
end

function [norm_L2,gradient_norm_L2] = fun_opt_persp_L2_only_h(I,pt_in_poly,I_pt,I_vect,x_param)
	if nargout > 1
		[err,J_err] = fun_opt_persp(I,pt_in_poly,I_pt,I_vect,x_param);
		[norm_L2,gradient_norm_L2] = error_norm_L2(err,J_err);
		gradient_norm_L2 = gradient_norm_L2(1,1:8);
	else
		[err] = fun_opt_persp(I,pt_in_poly,I_pt,I_vect,x_param);
		[norm_L2] = error_norm_L2(err);
	end
end


%  This version of the optimizer for the perspectivity enforce monotonic behavior by summing square of a_vect
%	to create r_vect spline control points
function [H_opt_normalized,r_opt] = optimize_perspectivity_monotonic(I,pt_in_poly,I_pt,I_vect,...
		invK,x_init,H_normalized,NB_R,R_vect)

	% We are using similar setting and enforce levenberg-marquardt
	options_perspective = optimoptions('lsqnonlin');
        options_perspective.Algorithm = 'levenberg-marquardt';
        options_perspective.Display = 'off';
        options_perspective.StepTolerance = 1e-12;
        options_perspective.FunctionTolerance = 1e-10;
        options_perspective.MaxIter = 200;
        options_perspective.SpecifyObjectiveGradient = false;
        options_perspective.CheckGradients = false;
        options_perspective.MaxFunctionEvaluations = 10000;

	% Sample the point use in the optimisation at least by 2
	permutation_0 = randperm(length(I_pt));
	% Select first half point
	keep = 0.3;
	n_half = floor(length(I_pt)*keep);
	rand_half = permutation_0(1:n_half);
	I_pt = I_pt(rand_half,:);
	pt_in_poly = pt_in_poly(rand_half,:);
	

	% Use sqrt of the additive component so that we do not need explicit constraint
	%	to keep monotonic spline
	%as_vect_init = sqrt(x_init(5:end));
	as_vect_init = sqrt(x_init(5:end));
	x_init = [x_init(1:4);as_vect_init];
	fun_opt = @(param) evaluate_error_homography_monotone_proxy_only_x_spline(I,pt_in_poly,I_pt,I_vect,...
		vector_to_perspectivity_minimal_param(param(1:4),invK),param(5:(4+NB_R)));

	% Calculate the difference between jacobian
	
	tic;
        [x_opt,~,~,~,output] = lsqnonlin(fun_opt,x_init,[],[],options_perspective);
        toc

	[I_pt_proj,r_proj,err] = evaluate_error_homography_monotone_detailed(I,pt_in_poly,I_pt,I_vect,...
        vector_to_perspectivity_minimal_param(x_opt(1:4),invK),x_opt(5:(4+NB_R)).^2);
        H_opt_normalized = vector_to_perspectivity_minimal_param_display(x_opt(1:4),invK)*H_normalized;
        r_opt = R_vect*x_opt(5:(4+NB_R)).^2;

end

function [H_opt_normalized,r_opt] = optimize_perspectivity(I,pt_in_poly,I_pt,I_vect,...
		invK,x_init,H_normalized,NB_R,R_vect,lb_angle,ub_angle)

	options_perspective = optimoptions('lsqnonlin');
	options_perspective.Algorithm = 'trust-region-reflective';
	options_perspective.Display = 'iter-detailed';
	options_perspective.StepTolerance = 1e-12;
	options_perspective.FunctionTolerance = 1e-10;
	options_perspective.MaxIter = 250;
	options_perspective.SpecifyObjectiveGradient = false;
	options_perspective.CheckGradients = false;
	options_perspective.MaxFunctionEvaluations = 20000;
	if nargin<12
		lb_angle = [-Inf;0];
		ub_angle = [Inf;pi/2];
	end
	%lb = [lb_angle;-Inf*ones(3,1);zeros(NB_R-1,1)];
	lb = [lb_angle;-Inf*ones(3,1);zeros(NB_R-1,1)];
	ub = [ub_angle;Inf*ones(2+NB_R,1)];
	fun_opt = @(param) evaluate_error_homography_monotone(I,pt_in_poly,I_pt,I_vect,...
	vector_to_perspectivity_minimal_param(param(1:4),invK),param(5:(4+NB_R)));
	% Calculate the initial point for the optimization

	%h_vect_init_not_normalized = homography_to_vector_perspective(...
	%vector_to_homography_perspective(h_vect_init)*inv(H_normalized));
	%[v_init_persp,scale] = vector_init_to_perspectivity_minimal_param(h_vect_init_not_normalized,invK);
	%x_init = [v_init_persp;a_vect_init./(scale.^2)];
	r_vect_init = R_vect*x_init(5:end);

	% To check that our initial point is correct, display the evalutation of the
	% photometric error
	%[I_pt_init,r_proj_init,err_init] = ...
	%        evaluate_error_homography_monotone_detailed(I,pt_in_poly,I_pt,I_vect,...
	%        vector_to_perspectivity_minimal_param(x_init(1:4),invK),x_init(5:end));
	%display_spline_function(r_vect_init,r_proj_init,I_vect,I_pt,...
	%                'Display spline interpolation initialization');
	tic;
	% When working with the constraints the value of the residual in optimization debug log
	%	is higher, this should not be a problem though
	[x_opt,~,~,~,output] = lsqnonlin(fun_opt,x_init,lb,ub,options_perspective);
	%x_opt = x_init;
	toc
	% Calculate the projected intensity we would want

	[I_pt_proj,r_proj,err] = evaluate_error_homography_monotone_detailed(I,pt_in_poly,I_pt,I_vect,...
	vector_to_perspectivity_minimal_param(x_opt(1:4),invK),x_opt(5:(4+NB_R)));
	H_opt_normalized = vector_to_perspectivity_minimal_param_display(x_opt(1:4),invK)*H_normalized;
	r_opt = R_vect*x_opt(5:(4+NB_R));
end

function [H_opt_normalized,r_opt,I_vect_opt] = optimize_perspectivity_full_cp(I,pt_in_poly,I_pt,I_vect,...
		invK,x_init_0,H_normalized,NB_R,R_vect,lb_angle,ub_angle)

	options_perspective = optimoptions('lsqnonlin');
	options_perspective.Algorithm = 'trust-region-reflective';
	options_perspective.Display = 'off';
	options_perspective.StepTolerance = 1e-12;
	options_perspective.FunctionTolerance = 1e-10;
	options_perspective.MaxIter = 250;
	options_perspective.SpecifyObjectiveGradient = false;
	options_perspective.CheckGradients = false;
	options_perspective.MaxFunctionEvaluations = 20000;
	if nargin<12
		lb_angle = [-Inf;0];
		ub_angle = [Inf;pi/2];
	end
	lb = [lb_angle;-Inf*ones(3,1);zeros(NB_R-1,1);zeros(NB_R,1)];
	ub = [ub_angle;Inf*ones(2+NB_R,1);255*ones(NB_R,1)];
	fun_opt = @(param) evaluate_error_homography_monotone(I,pt_in_poly,I_pt,transpose(param((5+NB_R):4+2*NB_R)),...
	vector_to_perspectivity_minimal_param(param(1:4),invK),param(5:(4+NB_R)));
	x_init = [x_init_0;transpose(I_vect)];
	tic;
	[x_opt,~,~,~,output] = lsqnonlin(fun_opt,x_init,lb,ub,options_perspective);
	toc
	% Calculate the projected intensity we would want
	[I_pt_proj,r_proj,err] = evaluate_error_homography_monotone_detailed(I,pt_in_poly,I_pt,I_vect,...
	vector_to_perspectivity_minimal_param(x_opt(1:4),invK),x_opt(5:(4+NB_R)));
	H_opt_normalized = vector_to_perspectivity_minimal_param_display(x_opt(1:4),invK)*H_normalized;
	r_opt = R_vect*x_opt(5:(4+NB_R));
	I_vect_opt = transpose(x_opt((5+NB_R):4+2*NB_R));
end

