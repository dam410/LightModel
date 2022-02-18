function [N1,N2] = vanishing_line_paper(Q)

	Q = Q*sign(det(Q))/(abs(det(Q))^(1/3));

	trq = trace(Q);
	triq = trace(inv(Q));

	a = (1/2)*(2*trq^3-9*trq*triq+27)/((trq^2-3*triq).^(3/2));
	%a = (1/2)*(trq^3-9*trq*triq*+27)/((trq^2-3*triq).^(3/2))
	alpha = 1/3*acos(a);

	l12 = sqrt(trq^2-3*triq)*(cos(alpha)-sqrt(3)/3*sin(alpha));
	l23 = sqrt(trq^2-3*triq)*2*sqrt(3)/3*sin(alpha);

	% Singular values are sorted in ascending order in Matlab
	[V,D] = eig(Q);
	V1 = V(:,1);
	V3 = V(:,3);
	N1 = sqrt(l12)*V1/norm(V1)+sqrt(l23)*V3/norm(V3);
	N1 = N1/norm(N1);
	N2 = sqrt(l12)*V1/norm(V1)-sqrt(l23)*V3/norm(V3);
	N2 = N2/norm(N2);
	
end
