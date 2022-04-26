% Create the data for real experiement
function [data,polys_2D,I_undistorted] = get_data_detect_real(img_file,cameraParams,tagSize,candleHeight,detection_type)
	% Initialize the data
	data = struct();
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
	K = transpose(cameraParams.Intrinsics.IntrinsicMatrix);
	[meshes_vector,S] = detect_groundtruth_scene(I_undistorted,K,...
		'CandleHeight',candleHeight,...
		'TagSize',tagSize);
	[polys_2D,polys_3D] = project_mesh(meshes_vector,T_cam,K);
	[polys] = mesh2polygons(meshes_vector,T_cam);


	[polys] = mesh2polygons(meshes_vector,T_cam);
	scene_plane_orientation = {};
	scene_plane_distance = {};
	scene_plane_distance_to_source = {};
	scene_plane_polygon = {};

	for i_poly=1:size(polys,1)
		i = i_poly;
		poly = polys(i_poly,:);
		scene_plane_orientation{i} = poly{1};
		scene_plane_distance{i} = poly{2};
		scene_plane_polygon{i} = polys_3D{i_poly};
		if (scene_plane_orientation{i}(3,3)<0)
			scene_plane_orientation{i} = scene_plane_orientation{i}*diag([1,-1,-1]);
			scene_plane_distance{i} = -poly{2};
		end
		scene_plane_distance_to_source{i} = ...
			-(transpose(scene_plane_orientation{i}(:,3))*S + scene_plane_distance{i});
	end
	
	% Create groundtruth data for scene elements pose parameters
	groundtruth = struct();
	groundtruth.SourcePosition = S;
	groundtruth.ScenePlaneOrientation = scene_plane_orientation;
	groundtruth.ScenePlanePosition = scene_plane_distance;
	groundtruth.ScenePlaneDistanceSource = scene_plane_distance_to_source;
	groundtruth.ScenePlanePolygon = scene_plane_polygon;
	data.groundtruth = groundtruth;
	data.K = K;
	data.T_cam = T_cam;

	% Warning this may be appeared to be redundancy with K, but here
	%	we also have the distorsion that we use to undistort the image
	data.Intrinsics = cameraParams;

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
