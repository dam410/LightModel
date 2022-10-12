function [data,I_rgb] = get_data_ex1_pure_synth(datafile,noise,nbIso)

	% Number of points calculated in the isocontours
	NB_SAMPLES = 10000;
	NB_KEPT = 2000;

	% Initialize the data
	data = struct();
	% Read the json file to extract data
	[mesh_,camera_,light_,render_] = decode_json_blender(datafile);
	% Import the camera pose and intrinsic
	[K,O,T_cam,h,l] = get_camera_matrices(camera_(1),render_(1));
	H_image = H_image_fcn(h,l);
	inv_T_cam = inv(T_cam);
	% Calculate the source
	S_4 = inv(T_cam)*light_(1).pose*[0;0;0;1];
	S = S_4(1:3)/S_4(4);
	D_2 = light_(1).data.distance.^2;
	coef_inv_quad = light_(1).data.data.quadratic_attenuation;
	[polys] = mesh2polygons(mesh_,T_cam);
	% Load the image in Matlab
	if nargin < 3
		nbIso = 2;
	end
	% Calculate isocontour on scene planes
	%	Calculate the polygon on image projection of the scene surface.
	K_im = [0,1,0;1,0,0;0,0,1]*H_image*K;
	[polys_2D,polys_3D] = project_mesh(mesh_,T_cam,K_im);
	%	Calculate two isocontours for each polygon surfaces
	isocontour_curve_params = {};
	isocontour_pts = {};
	scene_plane_orientation = {};
	scene_plane_distance = {};
	scene_plane_distance_to_source = {};
	scene_plane_polygon = {};
	% Consider all the planes on the scene
	for i=1:size(polys,1)
		poly = polys(i,:);
		scene_plane_orientation{i} = poly{1};
		scene_plane_distance{i} = poly{2};
		if (scene_plane_orientation{i}(3,3)<0)
			scene_plane_orientation{i} = scene_plane_orientation{i}*diag([1,-1,-1]);
			scene_plane_distance{i} = -poly{2};
		end
		scene_plane_distance_to_source{i} = ...
			-(transpose(scene_plane_orientation{i}(:,3))*S + scene_plane_distance{i});
		% In this file instead of detecting 2D points on isocontour,
		% we will calculate random points that belong to the real curve

		Y_source = transpose(polys{i,1})*S;
		X_source = polys{i,1}*[Y_source(1);Y_source(2);-polys{i,2}];
		inside_poly = inpolygon(X_source(1),X_source(2),...
			polys_3D{i}(1,:),polys_3D{i}(2,:));
		r_polygon = sqrt( (polys_3D{i}(1,:)-X_source(1)).^2 + ...
			(polys_3D{i}(2,:)-X_source(2)).^2 ) ;
		% Chosing the radius of circles :
		%	* Case 1 : The source projection is inside the polygons
		%	Then we choose betwen radius between [0.1,min(distance)-0.1]
		%	* Case 2 : The source projection is not inside the polygons
		%	Then we take circles of radius between r_min and r_max
		if inside_poly
			r_min = 0;
			r_max = min(r_polygon);
		else
			r_min = min(r_polygon);
			r_max = max(r_polygon);
		end
		r_iso = r_min:(r_max-r_min)/(nbIso+1):r_max;
		r_iso = r_iso(2:(end-1));
		pts = cell(1,nbIso);
		curve_params = cell(1,nbIso);
		for i_iso = 1:nbIso
			theta = 0:(2*pi/NB_SAMPLES):(2*pi);
			X_iso = r_iso(i_iso)*cos(theta)+Y_source(1);
			Y_iso = r_iso(i_iso)*sin(theta)+Y_source(2);
			Z_iso = -polys{i,2}*ones(size(X_iso));
			W_iso = polys{i,1}*[X_iso;Y_iso;Z_iso];
			X_im = K_im*W_iso;
			X_im = X_im./X_im(3,:);
			indices = 1:length(theta);
			circ_pts_inside = inpolygon(X_im(1,:),X_im(2,:),...
			polys_2D{i}(1,:),polys_2D{i}(2,:));
			indices_inside = indices(circ_pts_inside);
			indices_samples = datasample(indices_inside,NB_KEPT,'Replace',true);
			X_im_inside = X_im(:,indices_samples);
			data_points = transpose(X_im_inside([2,1],:));
			% Adding noise on the points calculated
			data_points = data_points+noise*randn(size(data_points));
			[ellipse_param,~] = ellipseFromPoints(data_points);
			curve_params{i_iso} = ellipse_param;
			pts{i_iso} = data_points;
		end
		isocontour_pts{i} = pts;
		isocontour_curve_params{i} = curve_params;
		scene_plane_polygon{i} = polys_3D{i};
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
	data.K = K_im;
	data.T_cam = T_cam;
end
