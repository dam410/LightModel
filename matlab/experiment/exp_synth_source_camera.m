% Experiments on the distance between camera and source with synthetic data 

% Read the txt files with names of the data in json format
folder_path = '../data/exp_synth_source_camera/';
str_file = fopen([folder_path,'filenames_source_cam.txt'],'r');
% For each values of the distance between camera and source
for i_dis = 1:20
	distance_cam_source = (i_dis-1)/20;
	for i_samp = 1:10
		line = fgetl(str_file);
		if i_dis == 20 && i_samp == 10
			filename = [folder_path,line];
			% Load the datafile
			[data,polys_2D,I,I_rgb,~] = get_data_no_detection(filename,1);
			% Detection of the isocontours with homography model
			data_td = isocontour_detection(data,I_rgb,polys_2D,'Mode','top-down');
			% Detection of the isocontours with rotation model (colocalised)
			data_td_circular = isocontour_detection(data,I_rgb,polys_2D,'Mode','top-down-circular');
			% Estimating the pose with global model
			[dssp_global,psp_global,ps_global] = case_8_np(data_td,norm(data.groundtruth.SourcePosition));
			% Estimating the pose with colocalised model
			[dssp_colocalised,psp_colocalised,ps_colocalised] = case_8_np(data_td_circular,norm(data.groundtruth.SourcePosition));
			% Measure the error and store it to be displayed
			[err_h_global,err_orient_global,err_d_global,err_s_global] = ...
				eval_pose_n(dssp_global,psp_global,ps_global,data_td);
			[err_h_colocalised,err_orient_colocalised,err_d_colocalised,err_s_colocalised] = ...
				eval_pose_n(dssp_colocalised,psp_colocalised,ps_colocalised,data_td_circular);
			end
	end
end

fclose(str_file);
