function [norm_L2,gradient_norm_L2] = error_norm_L2(err,J_err)
	residual_L2 = err.^2;
	norm_L2 = sum(residual_L2,1);
	if nargout > 1 & nargin > 1
		% Only perform this for non zero values in J_err
		%percentage_filled = nnz(J_err)/(prod(size(J_err)))*100
		diff_residual = 2*err;
		[i_J_err,j_J_err,v_J_err] = find(J_err);
		J_residual_L2 = sparse(i_J_err,j_J_err,diff_residual(i_J_err).*v_J_err,size(J_err,1),size(J_err,2));
		gradient_norm_L2 = sum(J_residual_L2,1);
	end
end
