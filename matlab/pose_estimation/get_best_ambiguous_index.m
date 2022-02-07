% This function help calculate the best values of a parameters
%	when we have several estimation of it which all provide two ambiguous solution.
%	and the index corresponding to the good values
function [V_best] = get_best_ambiguous_index(V,n_est,n_dim,index)
	V_best = zeros(n_dim,1);
	for i = 1:n_est
		i_start = 1+(i-1)*n_dim;
		i_end = i*n_dim;
		V_i = V(i_start:i_end,index(i));
		V_best = V_best+V_i/n_est;
	end
end
