function [dssp,psp,ps] = complete_plane_pose(data,dssp,psp,ps)
	n_p = length(psp);
	for i_p = 1:n_p
		points_3D = data.groundtruth.ScenePlanePolygon{i_p};
		points_3D = [points_3D,points_3D(:,1)]; 
		points_3D_bary = mean(points_3D,2);
		% Calculate the plane pose
		psp{i_p}{1}(4) = -transpose(psp{i_p}{1}(1:3))*points_3D_bary;
		dssp{i_p}{1} = -psp{i_p}{1}(4);
	end
end

