% This function help calculate the best values of a parameters
%	when we have several estimation of it which all provide two ambiguous solution.
% We compare the values provided the tools and give the list of index giving for each pair
%	the one closest to the consensus value of the parameters.
% Compute the average at the end in V_best and the index of the chosen value in index_best
function [V_best,index_best] = get_best_ambiguous(V,n_est,n_dim,compfcn)
	% We have nb_est > 1
	V_start = V(1:n_dim,1:2);
	V_best = zeros(n_dim,1);
	index_best = [];
	for i = [2:n_est,1]
		i_start = 1+(i-1)*n_dim;
		i_end = i*n_dim;
		V_i_1 = V(i_start:i_end,1);
		V_i_2 = V(i_start:i_end,2);
		[val_1,best_1] = min([compfcn(V_i_1,V_start(:,1)),compfcn(V_i_1,V_start(:,2))]);
		[val_2,best_2] = min([compfcn(V_i_2,V_start(:,1)),compfcn(V_i_2,V_start(:,2))]);
		[~,i_best] = min([val_1,val_2]);
		index_best = [index_best,i_best];
		V_best = V_best+V(i_start:i_end,i_best)/n_est;
		V_start = V(i_start:i_end,:);
	end
	index_best = index_best([n_est,1:(n_est-1)]);
end
