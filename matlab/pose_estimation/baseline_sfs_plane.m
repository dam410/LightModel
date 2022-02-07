% Function which use the image and its polygonal masks to estimate the pose based on an SfS algorithm
function [dssp,psp,ps,meshes] = baseline_sfs_plane(data,I_undistorted,polys_2D,SCALE)
	if nargin < 4
		SCALE = 4;
	end
	n_p = length(polys_2D);
	NB_R = 5;
	% For each scene plane, calculate the position of points inside polygons to generate a mask
	ind_pt_all = cell(n_p,1);
	for i_p = 1:n_p
		[pt_in_poly,I_pt,pt_barycenter,I_vect,ind_pt] = img_points_from_poly(I_undistorted,polys_2D{i_p},NB_R);
		ind_pt_all{i_p} = ind_pt;
	end
	% Initialize the scene elements poses output in correct format
	psp = cell(1,n_p);
	ps = {[0;0;0]};
	dssp = cell(1,n_p);
	% Calculate the global mask
	mask = zeros(size(I_undistorted,1),size(I_undistorted,2));
	mask(cell2mat(ind_pt_all)) = 1;
	% Downscale the mask and image accordingly (Not enough memory)
	mask_ds = imresize(mask,1/SCALE);
	mask_ds = logical(mask_ds);
	I_undistorted_grayscale = imresize(I_undistorted,1/SCALE);
	I_undistorted_grayscale = double(rgb2gray(I_undistorted_grayscale))/255;
	% Apply the SIRFS algorithm
	size(mask_ds)
	size(I_undistorted_grayscale)
	output = SIRFS(I_undistorted_grayscale, mask_ds, [], '');
	% Prepare the output mesh
	meshes = cell(1,n_p);
	% Now integrate the normal for each scene plane individually
	for i_p = 1:n_p
		normal_map = imresize(output.normal,SCALE);
		% Recreate a mask with only considered scene plane patch
		mask = zeros(size(I_undistorted,1),size(I_undistorted,2));
		mask(ind_pt_all{i_p}) = 1;
		mask = logical(mask);
		K = [0,1,0;1,0,0;0,0,1]*data.K;
		% Save the normal map, the camera matrix and the mask
		filename = '/home/dam/Documents/PostDoc_Damien/LightModel/data/controlled_sequences_SIRFS/temp_normal_map.mat';
		save(filename,'normal_map','K','mask');
		% Apply the normal integration algorithm
		python_folder = '/home/dam/Documents/PostDoc_Damien/LightModel/python/NormalIntegration';
		current_folder = pwd;
		command_normal_integration = ['cd ',python_folder,...
			'; /home/dam/anaconda3/envs/ni/bin/python methods/perspective_five_point_plane_fitting.py --path ',filename,...
			'; rm ',filename];
		system(command_normal_integration);
		ptCloud = pcread('/home/dam/Documents/PostDoc_Damien/LightModel/data/controlled_sequences_SIRFS/est_surface_perspective_five_point_plane_fitting.ply');
		ptCloud_cleaned = clean_point_cloud(ptCloud,1e-6);
		meshes{i_p} = ptCloud_cleaned.Location;
		P = [ptCloud_cleaned.Location,ones(ptCloud_cleaned.Count,1)];
		[U,S,V] = svd(transpose(P)*P);
		plane_pose = V(:,4)/norm(V(1:3,4));
		% Check if the normal is well oriented
		if plane_pose(3)<0
			plane_pose = -plane_pose;
		end
		psp{i_p} = {plane_pose};
		dssp{i_p} = {-plane_pose(4)};
	end
end
