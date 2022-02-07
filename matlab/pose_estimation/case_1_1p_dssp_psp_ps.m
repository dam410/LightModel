function [dssp,psp,ps] = case_1_1p_dssp_psp_ps(data,i_surface)

	if nargin < 2
		i_surface = 1;
	end

	% Known parameters
	n_p = 1;
	i_p = 1;
	dssp = cell(1,n_p);
	psp = cell(1,n_p);
	dssp{i_p} = {data.groundtruth.ScenePlaneDistanceSource{i_surface}};
	psp{i_p} = {[data.groundtruth.ScenePlaneOrientation{i_surface}(:,3);data.groundtruth.ScenePlanePosition{i_surface}]};
	ps = {data.groundtruth.SourcePosition};
end
