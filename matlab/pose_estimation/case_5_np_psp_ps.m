function [dssp,psp,ps] = case_5_1p_psp_ps(data)


	% Known parameters
        n_p = length(data.groundtruth.ScenePlaneDistanceSource);
        i_p = 1;
        psp = cell(1,n_p);
        for i_p = 1:n_p
                psp{i_p} = {[data.groundtruth.ScenePlaneOrientation{i_p}(:,3);data.groundtruth.ScenePlanePosition{i_p}]};
        end
        ps = {data.groundtruth.SourcePosition};

	% Estimation
	dssp = cell(1,n_p);
	for i_p = 1:n_p 
		dssp{i_p} = {distance_source_plane(psp{i_p}{1}(1:3),ps{1},psp{i_p}{1}(4))};
	end
end
