%% Load the image
%filename = '/home/dam/Documents/PostDoc_Damien/LightModel/data/exp_illustration/candle_planes.JPG';
%
%% Load the intrinsics
%load('../data/calibration/matlab_params_7.mat');
%
%[data,polys_2D,I_undis] = get_data_no_marker_real(filename,cameraParams,6);
%
%intrinsics = cameraParams.Intrinsics;
%K = [intrinsics.FocalLength(1),0,intrinsics.PrincipalPoint(1);...
%	0,intrinsics.FocalLength(2),intrinsics.PrincipalPoint(2);0,0,1];
%% Undistort the image and save it
%
%save('data_exp_illustration.mat');
%disp('Save data with polygon selection');
%pause;
load('data_exp_illustration.mat');

[data_td] = isocontour_detection(data,I_undis,polys_2D,'Mode','top-down','Display','on');
data_td.polys_2D = polys_2D;

[dssp,psp,ps] = case_8_np(data_td,1);
[dssp,psp,ps] = visualize_results_2(data_td,dssp,psp,ps);

save('data_detection_exp_illustration.mat');
%load('data_detection_exp_illustration.mat');


[pts_cell,I_pt_cell,I_vect_cell] = refinement_format(data_td,I_undis,polys_2D);
[dssp_bu_pc,psp_bu_pc,ps_bu_pc,err_bu_pc] = ...
	global_dense_isocontours_optimization(data_td,pts_cell,I_pt_cell,I_vect_cell,...
	dssp,psp,ps);

save('data_full_opti_detection_exp_illustration.mat');
%load('data_full_opti_detection_exp_illustration.mat');

visualize_results_2(data_td,dssp_bu_pc,psp_bu_pc,ps_bu_pc);
