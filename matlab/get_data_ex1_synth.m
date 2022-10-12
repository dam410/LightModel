function [data,I_rgb] = get_data_ex1_synth(datafile,noise,nbIso)
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
	img_name = render_(1).data.filepath;
	I_rgb = imread(img_name);
	I = double(rgb2gray(I_rgb));
	if nargin>1 && noise>0
		I = I+noise*randn(h,l);
		I = max(1,I);
		I = min(255,I);
		I_int = uint8(I);
		I_rgb = uint8(zeros(h,l,3));
		I_rgb(:,:,1) = I_int;
		I_rgb(:,:,2) = I_int;
		I_rgb(:,:,3) = I_int;
	end
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
		try 
			poly = polys(i,:);
			scene_plane_orientation{i} = poly{1};
			scene_plane_distance{i} = poly{2};
			if (scene_plane_orientation{i}(3,3)<0)
				scene_plane_orientation{i} = scene_plane_orientation{i}*diag([1,-1,-1]);
				scene_plane_distance{i} = -poly{2};
			end
			scene_plane_distance_to_source{i} = ...
				-(transpose(scene_plane_orientation{i}(:,3))*S + scene_plane_distance{i});
			try
				[pts,levels,curve_params] = detect_isocontour(I,...
					'Polygon',polys_2D{i},...
					'NbIsocontour',nbIso,...
					'Fitting','ellipse');
				isocontour_pts{i} = pts;
				isocontour_curve_params{i} = curve_params;
				scene_plane_polygon{i} = polys_3D{i};
			catch exception
				disp(['Problem of detection for image: ',...
					render_(1).data.filepath]);
			end
		catch exception
			disp('Not enough planes on the scene');
		end
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
