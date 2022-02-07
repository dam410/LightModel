function [H_opt_normalized,r_opt] = optimize_rotation_only(I,pt_in_poly,I_pt,I_vect,...
		invK,x_init_0,NB_R,R_vect)

	options_perspective = optimoptions('lsqnonlin');
	options_perspective.Algorithm = 'trust-region-reflective';
	options_perspective.Display = 'off';
	options_perspective.StepTolerance = 1e-12;
	options_perspective.FunctionTolerance = 1e-10;
	options_perspective.MaxIter = 250;
	options_perspective.SpecifyObjectiveGradient = false;
	options_perspective.CheckGradients = false;
	options_perspective.MaxFunctionEvaluations = 20000;
	lb = [-Inf;0;-Inf;zeros(NB_R-1,1)];
	ub = [Inf;pi/2;Inf*ones(NB_R,1)];
	fun_opt = @(param) evaluate_error_homography_monotone(I,pt_in_poly,I_pt,I_vect,...
		vector_to_rotation(param(1:2),invK),param(3:(2+NB_R)));
	x_init = x_init_0;

	[I_pt_proj,r_proj,err] = evaluate_error_homography_monotone_detailed(I,pt_in_poly,I_pt,I_vect,...
		vector_to_rotation(x_init(1:2),invK),x_init(3:(2+NB_R)));

	r_opt = R_vect*x_init(3:(2+NB_R));
	%display_spline_function(r_opt,r_proj,I_vect,I_pt,'Display the variance before rotation optimisation');
	%pause;


	%[output] = fun_opt(x_init);
	%figure('Name','Show the error at init point');
	%plot(1:length(output),output,'-b');
	%pause;
	tic;
	[x_opt,~,~,~,output] = lsqnonlin(fun_opt,x_init,lb,ub,options_perspective);
	toc
	% Calculate the projected intensity we would want
	[I_pt_proj,r_proj,err] = evaluate_error_homography_monotone_detailed(I,pt_in_poly,I_pt,I_vect,...
		vector_to_rotation(x_opt(1:2),invK),x_opt(3:(2+NB_R)));
	H_opt_normalized = vector_to_rotation(x_opt(1:2),invK);
	r_opt = R_vect*x_opt(3:(2+NB_R));
end



