function [dssp,psp,ps] = case_1_np_dssp_psp_ps(data)

	% Known parameters
	n_p = length(data.groundtruth.ScenePlaneDistanceSource);
	dssp = cell(1,n_p);
	psp = cell(1,n_p);
	for i_p = 1:n_p
		dssp{i_p} = {data.groundtruth.ScenePlaneDistanceSource{i_p}};
		psp{i_p} = {[data.groundtruth.ScenePlaneOrientation{i_p}(:,3);data.groundtruth.ScenePlanePosition{i_p}]};
	end
	ps = {data.groundtruth.SourcePosition};
end
