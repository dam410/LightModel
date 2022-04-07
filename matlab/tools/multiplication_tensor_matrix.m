function [T_out] = multiplication_tensor_matrix(T_in,M)
	[n,m,p] = size(T_in);
	[p,l] = size(M);
	T_out = reshape(reshape(T_in,[n*m,p])*M,[n,m,l]);
end
