% Function to extract on matlab struct data the configuration described in json format
function [data,polys_2D,I,I_rgb,I_test] = get_data_no_detection(datafile,noise)
	% Initialize the data
        data = struct();
	% Read the json file to extract data
	[mesh_,camera_,light_,render_] = decode_json_blender(datafile);
	% Import the camera pose and intrinsic
	[K,O,T_cam,h,l] = get_camera_matrices(camera_(1),render_(1));
	H_image = H_image_fcn(h,l);
	S_4 = inv(T_cam)*light_(1).pose*[0;0;0;1];
        S = S_4(1:3)/S_4(4);
	% Calculate the source
	[polys] = mesh2polygons(mesh_,T_cam);
	% Load the image in Matlab
	img_name = render_(1).data.filepath;
	I_rgb = imread(img_name);
	I = double(rgb2gray(I_rgb));
	if nargin>1 && noise>0
		size(I)
                I = I+noise*randn(h,l);
                I = max(1,I);
                I = min(255,I);
                I_int = uint8(I);
                I_rgb = uint8(zeros(h,l,3));
                I_rgb(:,:,1) = I_int;
                I_rgb(:,:,2) = I_int;
                I_rgb(:,:,3) = I_int;
        end
	new_name = [img_name(1:end-4),'_noised_',num2str(noise),img_name(end-3:end)];
	imwrite(I_rgb,new_name);
	data.img_name = new_name;

	% If I_synth is requested in output parameter generate the equivalent synthetic image
	if nargout > 4
		[I_test,~,~,~,P_camera] = render_shading_isocontour(l,h,...
			'Surface','Meshes',...
			'Meshes',mesh_,...
			'Scattering','Phong',...
			'Light',light_,...
			'Camera',camera_,...
			'Renderer',render_,...
			'UseImageTransformation',1);
		%I_synth = permute(I_test,[2 1 3]);
	end
	% Calculate isocontour on scene planes
	%	Calculate the polygon on image projection of the scene surface.
	K_im = [0,1,0;1,0,0;0,0,1]*H_image*K;
	%K_im = [0,1,0;1,0,0;0,0,1]*H_image*K;
	[polys_2D,polys_3D] = project_mesh(mesh_,T_cam,K_im);
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
		scene_plane_polygon{i} = polys_3D{i};
        end
        % Create groundtruth data for scene elements pose parameters
        groundtruth = struct();
        groundtruth.SourcePosition = S;
        groundtruth.ScenePlaneOrientation = scene_plane_orientation;
        groundtruth.ScenePlanePosition = scene_plane_distance;
        groundtruth.ScenePlaneDistanceSource = scene_plane_distance_to_source;
        groundtruth.ScenePlanePolygon = scene_plane_polygon;
        data.groundtruth = groundtruth;
        data.K = K_im;
        data.T_cam = T_cam;
end
