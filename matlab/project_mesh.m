% This function calculates the projection in image plane of the 3D points in the polygon
function [polys_2D,polys_3D] = project_mesh(mesh_all,T_cam,K_im)
	% Maybe we should vectorize but for now let's just make it quick and readable
	inv_T_cam = inv(T_cam);
	polys_2D = {};
	polys_3D = {};
	for j=1:length(mesh_all)
		mesh_ = mesh_all(j);
		T_mesh = mesh_.pose;
		for i = 1:length(mesh_.data.polygons)
			% Transform the points into
			P = inv_T_cam*T_mesh*[transpose(mesh_.data.vertices(mesh_.data.polygons(i).vertices_index+1,:));ones(1,length(mesh_.data.polygons(i).vertices_index))];
			P_3D = P(1:3,:)./P([4,4,4],:);
			P_2D = K_im*P_3D;
			poly_2D = P_2D(1:2,:)./P_2D([3,3],:);
			polys_2D{end+1} = poly_2D;
			polys_3D{end+1} = P_3D;
		end
	end
end
