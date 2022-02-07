% Load the calibration parameters
load('cameraParamsEndoscope4.mat','cameraParams');


filelist = '../data/endoscope_medical_2/list_images_quantitative.txt';
filename = '../data/endoscope_medical_2/2021-12-15_142911_IMG0019.jpg';
[data,polys_2D,I_undistorted] = get_data_detect_real_endo(filename,cameraParams,5.3,'no_detection');

% Get the groundtruth, the working windows (polys_2D) and the noised image (std-dev = 1)
I_rgb = imread(data.img_name);
I_gray = double(rgb2gray(I_rgb));

% Application de filtre
I_gray = medfilt2(I_gray,[5,5]);

% Detection with top-down method _td
data_td = isocontour_detection(data,I_gray,polys_2D,'Mode','top-down','N_circ',3,'Display','off');
disp('Detection with top-down approach done');
% Detection with circular cone constraint
data_co_td = isocontour_detection(data,I_gray,polys_2D,'Mode','top-down-circular','N_circ',3,'Display','off');
disp('Detection with top-down approach dedicated to colocalized model done');

figure('Name','Display the detection with top-down approach');
subplot(2,1,1);
imshow(I_rgb);
hold on;
for i_p = 1:length(data_td.isocontour.CurveParameters)
	for i_iso = 1:length(data_td.isocontour.CurveParameters{i_p}{1})
		pts = data_td.isocontour.Points{i_p}{1}{i_iso};
		plot(pts(:,2),pts(:,1),'+r');
		E = param2ellipse(data_td.isocontour.CurveParameters{i_p}{1}{i_iso});
		E_t = [0,1,0;1,0,0;0,0,1]*E*[0,1,0;1,0,0;0,0,1];
		displayEllipse(E_t);
	end
end
title('Top-down detection with perspectivity parametrisation');
subplot(2,1,2);
imshow(I_rgb);
hold on;
for i_p = 1:length(data_co_td.isocontour.CurveParameters)
	for i_iso = 1:length(data_co_td.isocontour.CurveParameters{i_p}{1})
		pts = data_co_td.isocontour.Points{i_p}{1}{i_iso};
		plot(pts(:,2),pts(:,1),'+r');
		E = param2ellipse(data_co_td.isocontour.CurveParameters{i_p}{1}{i_iso});
		E_t = [0,1,0;1,0,0;0,0,1]*E*[0,1,0;1,0,0;0,0,1];
		displayEllipse(E_t);
	end
end
title('Top-down detection with only rotation parametrisation');


% Calculate the first estimation based on close form solution
[dssp_td,psp_td,ps_td] = case_8_np(data_td,norm(data.groundtruth.SourcePosition));
% 	Colocalized model in the detection
[dssp_co_td,psp_co_td,ps_co_td] = case_8_np(data_co_td,norm(data.groundtruth.SourcePosition));
[dssp_td_co,psp_td_co,ps_td_co] = case_7_co_ps(data_td);
% Colocalized model in the detection and the resolution
[dssp_co_td_co,psp_co_td_co,ps_co_td_co] = case_7_co_ps(data_co_td);


% Show results for our estimation
visualize_results_2(data_td,dssp_td,psp_td,ps_td);
visualize_results_2(data_td,dssp_td_co,psp_td_co,ps_td_co);
visualize_results_2(data_co_td,dssp_co_td,psp_co_td,ps_co_td);
visualize_results_2(data_co_td,dssp_co_td_co,psp_co_td_co,ps_co_td_co);

[~,err_orient_td,~,err_s] = eval_pose_n(dssp_td,psp_td,ps_td,data_td);
[~,err_orient_td_co] = eval_pose_n(dssp_td_co,psp_td_co,ps_td_co,data_td);
[~,err_orient_co_td] = eval_pose_n(dssp_co_td,psp_co_td,ps_co_td,data_co_td);
[~,err_orient_co_td_co] = eval_pose_n(dssp_co_td_co,psp_co_td_co,ps_co_td_co,data_co_td);

disp('Error with complete model');
err_orient_td

disp('Error with complete model but td detection adapted');
err_orient_co_td

disp('Error with colocalised model');
err_orient_td_co

disp('Error with colocalised model and td detection adapted');
err_orient_co_td_co


