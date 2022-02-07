function [dssp,psp,ps] = case_8_1p(data,i_surface)

	% Select the surface in the data 
        if nargin < 2
                i_surface = 1;
        end

	% Known parameters

	% Estimation			
	i_p = 1;
	n_p = 1;
	psp = cell(1,n_p);
	dssp = cell(1,n_p);
	[Xc,N] = plane_orientation_from_circular_contours(data.K,data.T_cam,...
		data.isocontour.CurveParameters{i_surface});
	if size(N,2)==2
		P_s_1 = source_plane_from_circle_center_orientation(data.K,Xc(:,1),N(:,1));
		P_s_2 = source_plane_from_circle_center_orientation(data.K,Xc(:,2),N(:,2));
		dssp{i_p} = {0};
                psp{i_p} = {[N(:,1);0],[N(:,2);0]};
                ps = {P_s_1,P_s_2};

	else
		P_s = source_plane_from_circle_center_orientation(data.K,Xc,N);
		dssp{i_p} = {0};
		psp{i_p} = {[N;0]};
		ps = {P_s};
	end
end
