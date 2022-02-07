function [dssp,psp,ps] = case_4_1p_dssp(data,i_surface)

	% Select the surface in the data 
        if nargin < 2
                i_surface = 1;
        end

	% Known parameters
	dssp = {{data.groundtruth.ScenePlaneDistanceSource{i_surface}}};

	% Estimation
	i_p = 1;
	[Xc,N] = plane_orientation_from_circular_contours(data.K,data.T_cam,...
		data.isocontour.CurveParameters{i_surface});
	if size(N,2)==2
		L_s_1 = source_line_from_distance_source_plane(data.K,Xc(:,1),N(:,1),dssp{i_p}{1});
		L_s_2 = source_line_from_distance_source_plane(data.K,Xc(:,2),N(:,2),dssp{i_p}{1});
		psp = {{[N(:,1);0],[N(:,2);0]}};
		ps = {L_s_1,L_s_2};
	else
		L_s = source_line_from_distance_source_plane(data.K,Xc,N,dssp{i_p}{1});
		psp = {{[N;0]}};
		ps = {L_s};
	end
end
