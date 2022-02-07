function [X,err] = closest_point_lines(L_cell)
	% We calculate the point X that minimize
	% for each lines (d,m) : || m-cross(X,d) ||
	b = [];
	A = [];
	for i=1:length(L_cell)
		[d,m] = plucker_matrix_to_dm(L_cell{i});
		m = m/norm(d);
		d = d/norm(d);
		b = [b;m];
		A = [A;cross_antisym(d)];
	end
	X = A\b;
	if rank(A)<3 & false
		A
		b
		L_cell{i}
	end
	err = sqrt(mean((A*X-b).^2));
end
