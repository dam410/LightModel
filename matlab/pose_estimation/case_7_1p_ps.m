function [dssp,psp,ps] = case_7_1p_ps(data,i_surface)

	% Select the surface in the data 
        if nargin < 2
                i_surface = 1;
        end

	% Known parameters
	ps = {data.groundtruth.SourcePosition};
	% Estimation
	n_p = 1;
	i_p = 1;
	dssp = cell(1,n_p);
	psp = cell(1,n_p);
	[Xc,N] = plane_orientation_from_circular_contours(data.K,data.T_cam,...
		data.isocontour.CurveParameters{i_surface});
	if size(N,2)==2
		dssp{i_p} = {distance_source_plane_from_circle_center_source(data.K,Xc(:,1),N(:,1),ps{1}),...
				distance_source_plane_from_circle_center_source(data.K,Xc(:,2),N(:,2),ps{1})};
		d_1 = distance_plane_camera(N(:,1),ps{1},dssp{i_p}{1});
		d_2 = distance_plane_camera(N(:,2),ps{1},dssp{i_p}{2});
		psp{i_p} = {[N(:,1);d_1],[N(:,2);d_2]};
	else
		dssp{i_p} = {distance_source_plane_from_circle_center_source(data.K,Xc,N,ps{1})}; 
		d = distance_plane_camera(N,ps{1},dssp{i_p}{1});
		psp{i_p} = {[N;d]};
	end
end
