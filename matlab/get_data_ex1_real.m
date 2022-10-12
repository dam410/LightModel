% Create the data for real experiement
function [data,I_test] = get_data_ex1_real(img_file,cameraParams,S,tagSize,nbPlan,nbIso)
	% Initialize the data
	data = struct();
	I = imread(img_file);
	h = size(I,1);
	l = size(I,2);
	% Load camera parameters matrix
	intrinsics = cameraParams.Intrinsics;
	K = [intrinsics.FocalLength(1),0,intrinsics.PrincipalPoint(1);0,intrinsics.FocalLength(2),intrinsics.PrincipalPoint(2);0,0,1];
	% Undistort the image
	I = undistortImage(I,intrinsics,"OutputView","same");
	meshes_vector = detect_apriltag_m(I,K,tagSize,nbPlan,1.8);
	T_cam = eye(4);
	[polys_2D,polys_3D] = project_mesh(meshes_vector,T_cam,K);
	%K_matlab = [0,1,0;1,0,0;0,0,1]*K;

	% Detect the isocontours
	[pts,levels,curve_detected] = detect_isocontour_real_hand(I,polys_2D,nbIso);
		
	%figure('Name','Display the isocontour points and the fitted ellipse');
	%imshow(I);
	%hold on;
	%for i_poly = 1:length(polys_2D)
	%	for i = 1:length(curve_detected(i_poly,:))
	%		E = [0,1,0;1,0,0;0,0,1]*param2ellipse(curve_detected{i_poly,i})*[0,1,0;1,0,0;0,0,1];
	%		plot(pts{i_poly,i}(:,2),pts{i_poly,i}(:,1),'+r');
	%		displayEllipse(E);
	%	end
	%end
	
	if nargout>1
		% Show the generated image to check if Matlab understand the image well
		[I_test,~,~,~,P_camera] = render_shading_isocontour(l,h,...
			'Surface','Meshes',...
			'Meshes',meshes_vector,...
			'Scattering','Phong',...
			'LightParameters',[S;0.6],...
			'CameraIntrinsic',K,...
			'UseImageTransformation',0);
		I_test = permute(I_test,[2 1 3]);
	end

	% Calculate groundtruth scene plane parameters
	T_cam = eye(4);
	[polys] = mesh2polygons(meshes_vector,T_cam);
        isocontour_curve_params = {};
        isocontour_pts = {};
        scene_plane_orientation = {};
        scene_plane_distance = {};
        scene_plane_distance_to_source = {};
        scene_plane_polygon = {};
        % Consider all the planes on the scene
	% Sort the plane from right to left
	min_x = [];
	for i = 1:size(polys,1)
		min_x = [min_x,min(polys_2D{i}(1,:))];
	end
	
	[~,index_poly_sorted] = sort(min_x,'descend');
        for i=1:size(polys,1)
		i_poly = index_poly_sorted(i);
		poly = polys(i_poly,:);
		scene_plane_orientation{i} = poly{1};
		scene_plane_distance{i} = poly{2};
		isocontour_pts{i} = pts(i_poly,:);
		isocontour_curve_params{i} = curve_detected(i_poly,:);
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
        isocontour = struct();
        isocontour.Points = isocontour_pts;
        isocontour.CurveParameters = isocontour_curve_params;
        data.groundtruth = groundtruth;
        data.isocontour = isocontour;
        data.K = K;
        data.T_cam = T_cam;
end
