function [dssp,psp,ps] = case_5_1p_psp_ps(data,i_surface)

	% Select the surface in the data 
        if nargin < 2
                i_surface = 1;
        end

	% Known parameters
        n_p = 1;
        i_p = 1;
        psp = cell(1,n_p);
        psp{i_p} = {[data.groundtruth.ScenePlaneOrientation{i_surface}(:,3);data.groundtruth.ScenePlanePosition{i_surface}]};
        ps = {data.groundtruth.SourcePosition};

	% Estimation
        dssp = cell(1,n_p);
	dssp{i_p} = {distance_source_plane(psp{i_p}{1}(1:3),ps{1},psp{i_p}{1}(4))};
end
