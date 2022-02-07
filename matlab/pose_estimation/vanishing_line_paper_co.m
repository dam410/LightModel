% Using similar form as in the paper
% Specific to colocalised PLS and optical center
function [N] = vanishing_line_paper_co(Q)

	Q = Q*sign(det(Q))/(abs(det(Q))^(1/3));

	trq = trace(Q);
	triq = trace(inv(Q));

	% Singular values are sorted in ascending order in Matlab
	% Specific : smallest eigenvalue is the direction vector of the conic
	[V,D] = eig(Q);
	[~,index_sort] = sort(diag(D),'descend');
	V1 = V(:,index_sort(1));
	V3 = V(:,index_sort(3));

	% If colocalised a = -1
	a = (1/2)*(2*trq^3-9*trq*triq+27)/((trq^2-3*triq).^(3/2));
	%a = (1/2)*(trq^3-9*trq*triq*+27)/((trq^2-3*triq).^(3/2))
	if (1-norm(a))<1e-7
		%a = -1 ou 1
		l12 = 1/2*sqrt(trq^2-3*triq);
		N = sqrt(l12)*V1/norm(V1);
		N = N/norm(N);
	else
		alpha = 1/3*acos(a);
		l12 = sqrt(trq^2-3*triq)*(cos(alpha)-sqrt(3)/3*sin(alpha));
		l23 = 2/3*sqrt(3*trq^2-9*triq)*sin(alpha);
		N1 = sqrt(l12)*V1/norm(V1)+sqrt(l23)*V3/norm(V3);
		N1 = N1/norm(N1);
		N2 = sqrt(l12)*V1/norm(V1)-sqrt(l23)*V3/norm(V3);
		N2 = N2/norm(N2);
		N = [N1,N2];
	end
	
end
