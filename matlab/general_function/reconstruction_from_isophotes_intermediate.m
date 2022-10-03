% Function  that goes through all step from image to 3D reconstruction for one image
%	This function can generate all the intermediate results image or data
%	Parameters include
% Only work for real images, adding synthetic camera and images support later
%
function [] = reconstruction_from_isophotes(options)
	arguments
                options.imgName (1,1) string
		options.SaveName (1,1) string
		options.cameraParameters (1,1) cameraParameters
		% Parameters relative to surface patch extraction
		%	Number of patches considered in image, 0 if automatically detected with encoding
		options.NbPatches (1,1) double = 0
		%	Method of selection, 'auto' use encoding marker, 'hand' ask user
		%	'semi-auto' use encoding marker for GT but hand selection for
		options.MethodPatchSelection (1,1) string = 'auto'
		% 		Nature of light source (necessary for goundtruth)
		%			'candle' or 'endoscope' 
		options.LightSourceType (1,1) string
		%		Height of the Candle if auto mode
		options.CandleHeight (1,1) double
		%		Size of the Marker if auto mode
		options.TagSize (1,1) double
		%	Shrink applied to remove possible border of the patch to avoid isophote discontinuities
		options.Shrink (1,1) double = 3
		% Parameters relative to isophote extraction
		%	Number of isophote extracted by patches
		options.NbIsocontour (1,1) double = 2
		% 	Detection method used: 
		%		1) Bottom-up based on isolevel extraction in intensity field
		%		2) Top-Down based on dense optimisation using an empirical model on the
		%			intensity distribution
		options.DetectionMethod (1,1) string = 'top-down'
		%	Model used for the isophote projection, perspectivity, homography or
		%	pure rotation (endoscopic case)
		options.IsophoteModel (1,1) string = 'perspectivity'
		%	A cubic monotonic spline relates the radius to intensity level in image
		%	it is parametered by:
		%		Number of spline points used in model
		options.NbSplinePoint (1,1) double = 5;
		%	A model relates the coordinate in image to radius distance to the brightest point
		%	Gaussian blur parameters to smooth the image
                options.StandardDeviation (1,1) double = 3
                options.FilterSize (1,1) uint8  = 11
		% Parameters relative to first scene elements estimation
		%	Case studied what is a prior if groundtruth is used
		options.CaseNumber (1,1) double = 1
		% 	Endoscopic case
		options.EndoscopicCase (1,1) double
		%	Scale ?
		options.ReconstructionScale (1,1) double
		options.AutoScaleReconstruction (1,1) double = 0
		% Parameters relative to refinement ('photometric or geometric')
		options.RefinementCriteria (1,1) string = 'photometric'
		% Parameters relative to intermediate results
		options.ResultsFolder (1,1) string
        end
	
	resFolder = char(options.ResultsFolder);
	savename = char(options.SaveName);
	
	% Load the image 
	% Input:
	%	string options.img_name
	% Output:
	%	matrix n*m*3 I_rgb
	%	string img_file
	img_file = char(options.imgName);
	I_rgb = double(imread(options.imgName));

	% Pre-process the image 
	% Input:
	%	struct option.cameraParameters
	% Output:	 
	%	matrix n*m*3 I_rgb_undistorted
	%	matrix 3*3 K
	%	
	%  Undistorsion using intrinsics provided
	intrinsics = options.cameraParameters.Intrinsics;
	I_rgb_undistorted = undistortImage(I_rgb,intrinsics,"OutputView","same");
	end_img_file = length(img_file);
	new_name = [img_file(1:(end_img_file-4)),'_undistorted',img_file((end_img_file-3):end_img_file)];
	imwrite(I_rgb_undistorted/255,new_name);
	K = transpose(intrinsics.IntrinsicMatrix);
	save([resFolder,savename,'_data_image_undistorted.mat'],...
		'new_name','K','intrinsics');
	% Remove a vignetting effect

	% Extract patches in the image, automatically or manually by hand
	% Input:
	%	if 'auto' or 'semi-auto'
	%	double candleHeight
	%	double tagSize
	% Output:
	%	struct data groundtruth
	%	arrays double polys_2D	
	switch options.MethodPatchSelection
		case 'auto'
			if strcmp(options.LightSourceType,'candle')
				[meshes_vector,S] = detect_groundtruth_scene(I_rgb_undistorted,K,...
				'CandleHeight',options.CandleHeight,...
				'TagSize',options.TagSize);
			elseif strcmp(options.LightSourceType,'endoscope')
				[meshes_vector,S] = detect_groundtruth_scene_endo(I_rgb_undistorted,K,...
				'TagSize',options.TagSize,'Polygon',1);
				S = [0;0;0];
			end
			T_cam = eye(4);
			[polys_2D,polys_3D] = project_mesh(meshes_vector,T_cam,K);
			[polys] = mesh2polygons(meshes_vector,T_cam);	
			data = mesh2data(meshes_vector,S,T_cam,K,options.cameraParameters,new_name);
		otherwise
			disp('TODO: do something if no auto detection hand selection');
	end
	save([resFolder,savename,'_data_patch.mat'],'data','polys_2D');

	%imshow(I_rgb_undistorted);
	%hold on;
	%plot(polys_2D{1}(1,:),polys_2D{1}(2,:));
	%pause;
			
	% Extract isophotes from images on the given patches
	% Input:
	%	struct data
	%	matrix I_undistorted
	%	cell array polys_2D
	% Output:
	%	struct data groundtruth with detected isophote

	I_double = double(rgb2gray(uint8(I_rgb_undistorted)));
	switch options.DetectionMethod 
		case 'bottom-up'
			disp('TODO add bottom up part in the code');
		case 'top-down'
			[homog_p,~,spline_param] = dense_isocontours_detection(I_double,data.K,polys_2D,...
				'Mode',options.IsophoteModel,...
				'Display','on',...
				'Nb_R',options.NbSplinePoint);
			data = add_isocontours_to_homog(data,polys_2D,homog_p,2);
			save('all_data_temp.mat');
			%% Display the spline for each curve
			for i_p = 1:length(homog_p)
				% Get the spline parameters
				sparam = spline_param{i_p};
				r_opt = transpose(sparam(1,:));
				I_vect = sparam(2,:);
				H_opt = homog_p{i_p};
				R_vect = tril(ones(options.NbSplinePoint));
				inv_R_vect = inv(R_vect);
				r_opt = transpose(sparam(1,:));
				% Project the points
				[pt_in_poly,I_pt,pt_barycenter,I_vect,ind_pt] = img_points_from_poly(I_double,polys_2D{i_p},...
					options.NbSplinePoint);
				[I_pt_proj,r_proj,err_opt] = evaluate_error_homography_monotone_detailed(...
					I_double,pt_in_poly,I_pt,I_vect,H_opt,inv_R_vect*r_opt);
				% Display the spline
				[fig_handle] = display_spline_function(r_opt,r_proj,I_vect,I_pt,['Spline model fitted to patch n°',num2str(i_p)]);
				% Save the spline
				saveas(fig_handle,[resFolder,savename,'_figure_spline_fitted_patch_',num2str(i_p),'.png']);
			end
		otherwise
			disp('Unrecognized detection method');
	end
	save([resFolder,savename,'_data_detection.mat'],'data','homog_p','spline_param');

	% Display a figure with all the isophotes
	fig_conic = figure('Name','Fitted isophote');
	imshow(I_double/255);
	hold on;
	for i_p = 1:length(homog_p)
		displayEllipse(data.isocontour.CurveParameters{i_p}{1}{1},'g');
		displayEllipse(data.isocontour.CurveParameters{i_p}{1}{end},'g');
		plot(data.isocontour.Points{i_p}{1}{1}(:,2),data.isocontour.Points{i_p}{1}{1}(:,1),'+r');
		plot(data.isocontour.Points{i_p}{1}{end}(:,2),data.isocontour.Points{i_p}{1}{end}(:,1),'+r');
	end
	xlim([1,size(I_double,2)]);
	ylim([1,size(I_double,1)]);
	saveas(fig_conic,[resFolder,savename,'_figure_detected_isophotes.png']);

	% Estimate the scene elements parameters depending on the configuration
	% Input:
        %       struct data
        %       double scale (optional)
        %       cell array polys_2D
        % Output:
	%	cell array dssp
	%	cell array psp
	%	cell array ps
	switch options.CaseNumber
		case 1
			[dssp,psp,ps] = case_1_np_dssp_psp_ps(data)
		case 2
			[dssp,psp,ps] = case_2_np_dssp_psp(data)
		case 3
			[dssp,psp,ps] = case_3_np_dssp_ps(data)
		case 4
			[dssp,psp,ps] = case_4_np_dssp(data)
		case 5
			[dssp,psp,ps] = case_5_np_psp_ps(data)
		case 6
			[dssp,psp,ps] = case_6_np_psp(data)
		case 7
			if isfield(options,'EndoscopicCase')
				[dssp,psp,ps] = case_7_co_ps(data)
			else
				[dssp,psp,ps] = case_7_np_ps(data)
			end
		case 8
			if isfield(options,'ReconstructionScale')
				[dssp,psp,ps] = case_8_np_scaled(data,options.ReconstructionScale)
			elseif options.AutoScaleReconstruction
				[dssp,psp,ps] = case_8_np(data,norm(data.groundtruth.SourcePosition));
			else
				[dssp,psp,ps] = case_8_np(data)
			end
		otherwise
	end
	save([resFolder,savename,'_data_pose_estimation.mat'],'dssp','psp','ps');

	[~,~,~,fig_3D] = visualize_results_2(data,dssp,psp,ps);
	camup
	

	%% Possible global refinement of the scene elements parameters
	%% Input:
        %%       struct data
        %%       cell array dssp
        %%       cell array psp
        %%       cell array ps
	%%	matrix I_undistorted
	%%	cell array polys_2D
        %% Output:
        %%       cell array dssp_opt
        %%       cell array psp_opt
        %%       cell array ps_opt
	%%	double err_opt
	%[dssp,psp,ps] = complete_plane_pose(data,dssp,psp,ps);
	%if strcmp(options.RefinementCriteria,'photometric')
	%	[pts_cell,I_pt_cell,I_vect_cell] = refinement_format(data,I_double,polys_2D);
	%	[dssp_opt,psp_opt,ps_opt,err_opt] = global_dense_isocontours_optimization(...
	%		data,pts_cell,I_pt_cell,I_vect_cell,dssp,psp,ps);
	%else
	%	[dssp_opt,psp_opt,ps_opt,err_opt] = refine_pose_isocontours(...
	%		data,dssp,psp,ps,1);
	%end
	%save([resFolder,'reconstruction_from_isophotes_pose_estimation_opt.mat'],...
	%	'dssp_opt','psp_opt','ps_opt','err_opt');

	%visualize_results_2(data,dssp_opt,psp_opt,ps_opt);
	
end
