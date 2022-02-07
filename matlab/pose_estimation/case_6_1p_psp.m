function [dssp,psp,ps] = case_6_1p_psp(data,i_surface)

        % Select the surface in the data 
        if nargin < 2
                i_surface = 1;
        end

        % Known parameters
        n_p = 1;
        i_p = 1;
        psp = cell(1,n_p);
        psp{i_p} = {[data.groundtruth.ScenePlaneOrientation{i_surface}(:,3);data.groundtruth.ScenePlanePosition{i_surface}]};

	% Estimation
	[Xc,N] = plane_orientation_from_circular_contours(data.K,data.T_cam,...
		data.isocontour.CurveParameters{i_surface});
	if size(N,2)==2
		L_s_1 = source_line_from_plane_pose(data.K,Xc(:,1),psp{i_p}{1}(1:3),psp{i_p}{1}(4));
		L_s_2 = source_line_from_plane_pose(data.K,Xc(:,2),psp{i_p}{1}(1:3),psp{i_p}{1}(4));
		dssp = {{0}};
		ps = {L_s_1,L_s_2};
	else
		[L_s] = source_line_from_plane_pose(data.K,Xc,psp{i_p}{1}(1:3),psp{i_p}{1}(4));
		dssp = {{0}};
		ps = {L_s};
	end
	
end
