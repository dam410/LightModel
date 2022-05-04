load('../data/calibration/matlab_params_specular_raw.mat');
img_file = '/home/dam/Documents/PostDoc_Damien/LightModel/data/exp_specular_0/DSC_1312.png';
tagSize = 30;
candleHeight = 25;
detection_type = 'None';
[data,polys_2D,I_undistorted] = get_data_detect_real(img_file,cameraParams,tagSize,candleHeight,detection_type);

% Visualize the camera, the scene plane and the light source
[dssp,psp,ps] = case_1_np_dssp_psp_ps(data);
visualize_results_2(data,dssp,psp,ps);

% Simulation with Matlab

[I_test,I_raw,D_source,Lj,P_camera] = render_shading_isocontour(4748,3153,...
                'Surface','Plane',...
                'LightType','PLS',...
                'Scattering','Phong',...
		'ScatteringParameters',[0.9,0.4,550],...
                'SurfaceParameters',-psp{1}{1},...
                'LightParameters',[data.groundtruth.SourcePosition;0.6],...
                'CameraIntrinsic',data.K,...
		'UseImageTransformation',0);

figure('Name','Simulate image');
imshow(transpose(I_test(:,:,1)));



% Move to a new coordinate system centered on the plane with z=0
R_plane2cam = data.groundtruth.ScenePlaneOrientation{1};
P_markers_cam = data.groundtruth.ScenePlanePolygon{1};


% Position of the markers in scene plane frame (somme mm difference)
P_markers_plane = transpose(R_plane2cam)*P_markers_cam - transpose(R_plane2cam)*P_markers_cam(:,1);
% Position of the camera
C = -transpose(R_plane2cam)*P_markers_cam(:,1);
% Rotation to get into camera frame
R = R_plane2cam
% Position of the source in scene plane frame
S = transpose(R_plane2cam)*data.groundtruth.SourcePosition + C;
% Projection matrix into image plane
P = data.K*[R,C];



data_json = struct();
data_json.Source = transpose(S);
data_json.Markers = transpose(P_markers_plane);
data_json.cameraIntrinsicK = data.K;
data_json.RadialDistortion = cameraParams.Intrinsics.RadialDistortion;
data_json.TangentialDistortion = cameraParams.Intrinsics. TangentialDistortion;
data_json.CameraPosition = transpose(C);
data_json.RotationToCameraFrame = R;
data_json.ProjectionMatrix = P;
% Put data in JSON
fileID = fopen('specular_scene_0.json','w');
fprintf(fileID,jsonencode(data_json));
fclose(fileID);

