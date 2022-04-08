% Similar to the matlab function ppval, we evaluate a piecewise polynomial
%	function and calculate its jacobian in the point and the polynomial coefficients
% Input Parameters :
%	* pp is struct composed of the piecewise polynomial properties
%		'coefs' ->  (n_pieces * n_deg) matrix : polynomial coefficients sort by decreasing order
%		'breaks' -> (n_pieces+1) array : x value separator for the piecewise definition
%		'order' ->  degree-1 of the polynomial
%		'dim' -> 1 dimension
%		'form' -> 'pp'
%	* xx is the x valued data
function [v,J_v] = ppval_derivative(pp,J_pp_mat,J_breaks,xx,J_xx)
	n_pieces = pp.pieces;
	% Check for each xx values the interval where it belongs, as the function
	%	is supposed C1, this step is not considered in the Jacobian consideration
	[ ~,index] = histc(xx,[-inf,pp.breaks(2:n_pieces),inf]);
	% Centered the values in their interval
	xx_centered = xx-transpose(pp.breaks(index));
	J_xx_centered = J_xx-J_breaks(index,:);
	% Reshape the polynomial coefficients with the correct output piecewise
	coefs = pp.coefs(index,:);
	coefs_derivate = [0,(pp.order-1):-1:1].*pp.coefs(index,[1,1:(end-1)]);
	% Nested multiplication for the evaluation of the 
	v = coefs(:,1);
	u = coefs_derivate(:,1);
	J_a = permute(J_pp_mat(index,1,:),[1,3,2]);
	for i=2:pp.order
		v = xx_centered(:).*v + coefs(:,i);
		u = xx_centered(:).*u + coefs_derivate(:,i);
		J_a = diag(xx_centered)*J_a + permute(J_pp_mat(index,i,:),[1,3,2]);
	end
	J_v = diag(u)*J_xx_centered+J_a;
end
