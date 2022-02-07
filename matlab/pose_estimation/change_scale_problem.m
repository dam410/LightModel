function [dssp,psp,ps] = change_scale_problem(dssp,psp,ps,scale)
	n_p = length(dssp);
	ps{1} = scale*ps{1};
	% Then once the source is known, every other parameters can be calculated
	for i_p = 1:n_p
		dssp{i_p}{1} = scale*dssp{i_p}{1};
		psp{i_p}{1}(4) = scale*psp{i_p}{1}(4);
	end

end
