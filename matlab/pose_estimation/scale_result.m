function [dssp_out,psp_out,ps_out] = scale_result(dssp,psp,ps,norm_ps)
	nb_plane = length(psp);
	lambda = norm_ps/norm(ps{1});
	dssp_out = dssp;
	psp_out = psp;
	ps_out = ps;
	ps_out{1} = ps{1}*lambda;
	for i = 1:nb_plane
		dssp_out{i}{1} = lambda*dssp{i}{1};
		psp_out{i}{1}(4) = lambda*psp{i}{1}(4);
	end
end
