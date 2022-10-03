
% Load the image, undistort etc ....
load('cameraParamsEndoscope4.mat','cameraParams');

img_file = '/home/dam/Documents/PostDoc_Damien/LightModel/data/endoscope_medical_2/2021-12-15_142911_IMG0019.jpg';%2021-12-15_142824_IMG0017.jpg';

tagSize = 0.008;


% Load camera parameters matrix
intrinsics = cameraParams.Intrinsics;
K = [intrinsics.FocalLength(1),0,intrinsics.PrincipalPoint(1);...
0,intrinsics.FocalLength(2),intrinsics.PrincipalPoint(2);0,0,1];
T_cam = eye(4);
% Undistort the image
I = imread(img_file);
I_undistorted = undistortImage(I,intrinsics,"OutputView","same");
new_name = [img_file(1:end-4),'_undistorted',img_file(end-3:end)]
imwrite(I_undistorted,new_name);
data = struct();
data.img_name = new_name;

K = transpose(cameraParams.Intrinsics.IntrinsicMatrix);

filepath = 'temp.pgm';
imwrite(I_undistorted,filepath);
K = K*diag([1,1,1]);
[cell_R,cell_t,cell_id,cell_err] = detect_apriltag(filepath,K,tagSize);
% Find the scene elements using their id
meshes_vector = [];
S = [0;0;0];
[meshes_vector] = detect_groundtruth_scene_endo(I_undistorted,K,'Polygon',1);
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
data.Intrinsics = cameraParams;




% Now for each marker, we display in green iteratively, their number and we create polygonal roi
% that we display and save
figure('Name','Selection the polygon for markers');
imshow(I);

[dssp,psp,ps] = case_1_np_dssp_psp_ps(data);
[dssp,psp,ps] = visualize_results_2(data,dssp,psp,ps);
%
%
%%Sticker brown 3 -> cell_id = [20 21 10 25]
%[Lia,Locb] = ismember([26 28 29 4],cell_id);
%%[Lia,Locb] = ismember([0 5 27 16],cell_id);
%%%% Sticker black 1 -> cell_id = [0 5 27 16]
%%%[Lia,Locb] = ismember([0 5 27 16],cell_id);
%
%Locb = Locb(find(Locb>0));
%if length(find(Locb>0))>0
%	nb_marqueur = length(Locb);
%	for i=Locb
%		id_marker = cell_id(i);
%		P_marker = K*cell_t{i};
%		P_marker = P_marker/P_marker(3);
%		text(P_marker(1),P_marker(2),num2str(id_marker),'Color','g');
%	end
%	roi_stick_1 = drawpolygon;
%	roi_stick = roi_stick_1;
%	nb_pt = length(roi_stick.Position);
%	cell_R_i = cell_R(Locb);
%	cell_t_i = cell_t(Locb);
%	% Calculate the homography on the scene plane
%	%[mesh_poly] = marker2polygon(cell_R(Locb),cell_t(Locb),tagSize,1);
%	N_marqueurs = [];
%	P = [];
%	P_marqueur = [];
%	for j=1:nb_marqueur
%                R = transpose(cell_R_i{j});
%                t = cell_t_i{j};
%                N_marqueurs = [N_marqueurs,-R(:,3)];
%                P = [P,2*tagSize*R(:,1)+t,2*tagSize*R(:,2)+t,t];
%                P_marqueur = [P_marqueur,t];
%        end
%        Nd = estimate_plane(P);
%	R_plane = [null(Nd(1:3)*transpose(Nd(1:3))),Nd(1:3)];
%        if det(R_plane)<0
%                R_plane = R_plane(:,[2,1,3]);
%        end
%	% Estimate the coordinate on the plane of the polygon
%	P_poly_unscale = inv(K)*[transpose(roi_stick.Position);ones(1,nb_pt)];
%	poly_scale = -Nd(4)./(transpose(Nd(1:3))*P_poly_unscale);
%	P_poly = poly_scale.*P_poly_unscale;
%	% Calculate the coordinates of the points on the marker plane
%	p_poly = transpose(R_plane)*P_poly;
%	p_marqueur = transpose(R_plane)*P_marqueur;
%	% Calculate polypoints coordinates using marker points as base
%	% Use the first marker as origin
%	p_poly_base = p_poly(1:2,:) - p_marqueur(1:2,3);
%	p_marqueur_base = p_marqueur(1:2,1:2) - p_marqueur(1:2,3);
%	p_poly_coord = inv(p_marqueur_base)*p_poly_base;
%	% Get the mesh
%	[mesh_] = marker2polygon(cell_R_i,cell_t_i,tagSize,1);
%	% Recalculate the mesh_ with polygon points
%	mesh_.data.normals = transpose(p_poly_coord)*mesh_.data.normals(1:2,:);
%	mesh_.data.vertices = transpose(p_poly_coord)*(mesh_.data.vertices(1:2,:)-mesh_.data.vertices(3,:))+mesh_.data.vertices(3,:);
%	mesh_.data.polygons.vertices_index = 0:1:(length(p_poly_coord)-1);
%end
%% Find a way to visualize, the resulting new mesh
%[polys_2D,polys_3D] = project_mesh(mesh_,T_cam,K);
%[polys] = mesh2polygons([mesh_],T_cam);
%scene_plane_orientation = {};
%scene_plane_distance = {};
%scene_plane_distance_to_source = {};
%scene_plane_polygon = {};
%
%for i_poly=1:size(polys,1)
%	i = i_poly;
%	poly = polys(i_poly,:);
%	scene_plane_orientation{i} = poly{1};
%	scene_plane_distance{i} = poly{2};
%	scene_plane_polygon{i} = polys_3D{i_poly};
%	if (scene_plane_orientation{i}(3,3)<0)
%		scene_plane_orientation{i} = scene_plane_orientation{i}*diag([1,-1,-1]);
%		scene_plane_distance{i} = -poly{2};
%	end
%	scene_plane_distance_to_source{i} = ...
%		-(transpose(scene_plane_orientation{i}(:,3))*S + scene_plane_distance{i});
%end
%
%% Create groundtruth data for scene elements pose parameters
%data = struct();
%data.img_name = img_file;
%groundtruth = struct();
%groundtruth.SourcePosition = zeros(3,1);
%groundtruth.ScenePlaneOrientation = scene_plane_orientation;
%groundtruth.ScenePlanePosition = scene_plane_distance;
%groundtruth.ScenePlaneDistanceSource = scene_plane_distance_to_source;
%groundtruth.ScenePlanePolygon = scene_plane_polygon;
%data.groundtruth = groundtruth;
%data.K = K;
%data.T_cam = T_cam;
%data.Intrinsics = cameraParams;
%% Pose estimation
%[dssp,psp,ps] = case_1_np_dssp_psp_ps(data);
%[dssp,psp,ps] = visualize_results_2(data,dssp,psp,ps);
%%
%%
%%
%%%
%%%% Sticker black 2 -> cell_id = [26 28 29 4]
%%%[Lia,Locb] = ismember([26 28 29 4],cell_id);
%%%Locb = Locb(find(Locb>0));
%%%if length(find(Locb>0))>0
%%%	%meshes_vector = [meshes_vector,marker2polygon(cell_R(Locb),cell_t(Locb),tagSize,1.7)];
%%%        for i=Locb
%%%                id_marker = cell_id(i);
%%%                P_marker = K*cell_t{i};
%%%                P_marker = P_marker/P_marker(3);
%%%                text(P_marker(1),P_marker(2),num2str(id_marker),'Color','g');
%%%        end
%%%	roi_stick_2 = drawpolygon;
%%%end
%%%
%%%% Sticker brown 3 -> cell_id = [20 21 10 25]
%%%[Lia,Locb] = ismember([20 21 10 25],cell_id);
%%%Locb = Locb(find(Locb>0));
%%%if length(find(Locb>0))>0
%%%	%meshes_vector = [meshes_vector,marker2polygon(cell_R(Locb),cell_t(Locb),tagSize,1.7)];
%%%        for i=Locb
%%%                id_marker = cell_id(i);
%%%                P_marker = K*cell_t{i};
%%%                P_marker = P_marker/P_marker(3);
%%%                text(P_marker(1),P_marker(2),num2str(id_marker),'Color','g');
%%%        end
%%%	roi_stick_3 = drawpolygon;
%%%end
%%%
%%%%[meshes_vector,S] = detect_groundtruth_scene_endo(I_undistorted,K,options)
%%%
%%%
%%%function [poly] = 
%%%end
