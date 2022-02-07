function [d,m] = plucker_matrix_to_dm(L)
	d = -L(1:3,4);
	m = [-L(2,3);L(1,3);-L(1,2)];
end
