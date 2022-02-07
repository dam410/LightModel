function [norm_L1,gradient_norm_L1] = error_norm_L1(err,J_err)
	EPS_PSEUDO = 1e-4;	
	residual_L1 = sqrt(err.^2 + EPS_PSEUDO);
	norm_L1 = sum(residual_L1,1);
	if nargout > 1 & nargin > 1
		% Only perform this for non zero values in J_err
		%percentage_filled = nnz(J_err)/(prod(size(J_err)))*100
		diff_residual = err./residual_L1;
		[i_J_err,j_J_err,v_J_err] = find(J_err);
		J_residual_L1 = sparse(i_J_err,j_J_err,diff_residual(i_J_err).*v_J_err,size(J_err,1),size(J_err,2));
		%J_residual_L1 = diff_residual.*J_err;
		gradient_norm_L1 = sum(J_residual_L1,1);
	end
end
