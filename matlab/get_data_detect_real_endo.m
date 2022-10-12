% Create the data for real experiement
function [data,polys_2D,I_undistorted] = get_data_detect_real_endo(img_file,cameraParams,tagSize,detection_type)
	% Initialize the data
	I = imread(img_file);
	h = size(I,1);
	l = size(I,2);
	T_cam = eye(4);

	% Load camera parameters matrix
	intrinsics = cameraParams.Intrinsics;
	K = [intrinsics.FocalLength(1),0,intrinsics.PrincipalPoint(1);...
	0,intrinsics.FocalLength(2),intrinsics.PrincipalPoint(2);0,0,1];
	% Undistort the image
	I_undistorted = undistortImage(I,intrinsics,"OutputView","same");
        new_name = [img_file(1:end-4),'_undistorted',img_file(end-3:end)]
        imwrite(I_undistorted,new_name);

	K = transpose(cameraParams.Intrinsics.IntrinsicMatrix);
	[meshes_vector,S] = detect_groundtruth_scene_endo(I_undistorted,K,...
		'TagSize',tagSize);

	data = mesh2data(meshes_vector,S,T_cam,K,intrinsics,new_name);

	% How the isocontours are detected
	switch detection_type
		case 'perspectivity_photometric_optimization'
			I_double = double(rgb2gray(I_undistorted));
			[homog_p] = dense_isocontours_detection(I_double,data.K,polys_2D,...
				'Mode','perspectivity',...
				'Display','off',...
				'Nb_R',5);
			data = add_isocontours_to_homog(data,polys_2D,homog_p,2);
		case 'ellipse_detection'
			[pts,levels,curve_detected] = detect_isocontour_real(I_undistorted,polys_2D,4);
			isocontour_curve_params = {};
			isocontour_pts = {};
			for i=1:size(polys,1)
				isocontour_pts{i} = pts(i,:);
				isocontour_curve_params{i} = curve_detected(i,:);
			end
			isocontour.Points = isocontour_pts;
			isocontour.CurveParameters = isocontour_curve_params;
			data.isocontour = isocontour;
		case 'perspectivity_photometric_optimization_ambig'
			I_double = double(rgb2gray(I_undistorted));
			[homog_p_1,homog_p_2] = dense_isocontours_detection(I_double,...
				data.K,polys_2D,...
				'Mode','perspectivity',...
				'Display','off',...
				'Nb_R',5);
			data = add_isocontours_to_homogs_ambig(data,polys_2D,homog_p_1,homog_p_2,2);
		otherwise
			disp('Isocontours have not been detected');
	end
end
