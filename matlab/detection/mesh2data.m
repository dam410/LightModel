function [data] = mesh2data(meshes_vector,S,T_cam,K,cameraParams,img_name)
	
	data = struct();
	data.img_name = img_name;
	intrinsics = cameraParams.Intrinsics;

	[polys_2D,polys_3D] = project_mesh(meshes_vector,T_cam,K);
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
end
