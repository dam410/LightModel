% Create the data for real experiement
function [data,polys_2D,I] = get_data_no_detection_real(img_file,cameraParams,nbPlan)
	% Initialize the data
	data = struct();
	I = imread(img_file);
	h = size(I,1);
	l = size(I,2);
	% Load camera parameters matrix
	intrinsics = cameraParams.Intrinsics;
	K = [intrinsics.FocalLength(1),0,intrinsics.PrincipalPoint(1);...
		0,intrinsics.FocalLength(2),intrinsics.PrincipalPoint(2);0,0,1];
	% Undistort the image and save it
	I = undistortImage(I,intrinsics,"OutputView","same");
	new_name = [img_file(1:end-4),'_undistorted',img_file(end-3:end)]
	imwrite(I,new_name);
	data.img_name = new_name;
	% Ask user to draw nbPlan ROI
	polys_2D = cell(1,nbPlan);
	figure('Name','Image');
	imshow(I);
	for i=1:nbPlan
		poly = drawpolygon();
		polys_2D{i} = transpose(poly.Position);
	end
	%meshes_vector = detect_apriltag_m(I,K,tagSize,nbPlan,1.8);
	%T_cam = eye(4);
	%[polys_2D,polys_3D] = project_mesh(meshes_vector,T_cam,K);
	I = double(rgb2gray(I));

        %[polys] = mesh2polygons(meshes_vector,T_cam);
        %scene_plane_orientation = {};
        %scene_plane_distance = {};
        %scene_plane_distance_to_source = {};
        %scene_plane_polygon = {};
	%% Consider all the planes on the scene
        %% Sort the plane from right to left
        %min_x = [];
        %for i = 1:size(polys,1)
        %        min_x = [min_x,min(polys_2D{i}(1,:))];
        %end
        %[~,index_poly_sorted] = sort(min_x,'descend');
        %for i=1:size(polys,1)
        %        i_poly = index_poly_sorted(i);
        %        poly = polys(i_poly,:);
        %        scene_plane_orientation{i} = poly{1};
        %        scene_plane_distance{i} = poly{2};
        %        scene_plane_polygon{i} = polys_3D{i_poly};
        %        if (scene_plane_orientation{i}(3,3)<0)
        %                scene_plane_orientation{i} = scene_plane_orientation{i}*diag([1,-1,-1]);
        %                scene_plane_distance{i} = -poly{2};
        %        end
        %        scene_plane_distance_to_source{i} = ...
        %                -(transpose(scene_plane_orientation{i}(:,3))*S + scene_plane_distance{i});
        %end
	%polys_2D = polys_2D(index_poly_sorted);

        % Create groundtruth data for scene elements pose parameters
        %groundtruth = struct();
        %groundtruth.SourcePosition = S;
        %groundtruth.ScenePlaneOrientation = scene_plane_orientation;
        %groundtruth.ScenePlanePosition = scene_plane_distance;
        %groundtruth.ScenePlaneDistanceSource = scene_plane_distance_to_source;
        %groundtruth.ScenePlanePolygon = scene_plane_polygon;
        %data.groundtruth = groundtruth;
        data.K = K;
        data.T_cam = eye(4);
	% Warning this may be appeared to be redundancy with K, but here
	%	we also have the distorsion that we use to undistort the image
	data.Intrinsics = cameraParams;
end
