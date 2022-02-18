% Experiments on the distance between camera and source with synthetic data 

% Read the txt files with names of the data in json format
folder_path = '../data/exp_synth_source_camera/';
str_file = fopen([folder_path,'filenames_source_cam.txt'],'r');
% Size is 20 value for distance, 6 scene plane per scene multiply by 5 samples for scene gen
err_s_global_all = zeros(20,5);
err_orient_global_all = zeros(20,6*5);
err_orient_colocalised_all = zeros(20,6*5);
% For each values of the distance between camera and source
for i_dis = 1:20
	distance_cam_source = (i_dis-1)/20;
	for i_samp = 1:5
		line = fgetl(str_file);
		filename = [folder_path,line];
		% Load the datafile
		[data,polys_2D,I,I_rgb] = get_data_no_detection(filename,1);
		% Detection of the isocontours with homography model
		data_td = isocontour_detection(data,I_rgb,polys_2D,'Mode','top-down');
		% Detection of the isocontours with rotation model (colocalised)
		data_td_circular = isocontour_detection(data,I_rgb,polys_2D,'Mode','top-down-circular');
		% Estimating the pose with global model
		[dssp_global,psp_global,ps_global] = case_8_np(data_td,norm(data.groundtruth.SourcePosition));
		% Estimating the pose with colocalised model
		[dssp_colocalised,psp_colocalised,ps_colocalised] = case_7_co_ps(data_td_circular);
		% Measure the error and store it to be displayed
		[err_h_global,err_orient_global,err_d_global,err_s_global] = ...
			eval_pose_n(dssp_global,psp_global,ps_global,data_td);
		[err_h_colocalised,err_orient_colocalised,err_d_colocalised,err_s_colocalised] = ...
			eval_pose_n(dssp_colocalised,psp_colocalised,ps_colocalised,data_td_circular);
		% Visualize results for debugging only
		%visualize_results_2(data_td,dssp_global,psp_global,ps_global);
		%visualize_results_2(data_td_circular,dssp_colocalised,psp_colocalised,ps_colocalised);
		err_s_global_all(i_dis,i_samp) = err_s_global;
		err_orient_global_all(i_dis,(i_samp-1)*6+(1:6)) = err_orient_global;
		err_orient_colocalised_all(i_dis,(i_samp-1)*6+(1:6)) = err_orient_colocalised;
	end
end

fclose(str_file);
