function [L,err] = source_line_from_planes(planes)
	Ns = cell2mat(planes);
	M = Ns(1:3,:);
	A = M*transpose(M);
	try
		[U,S,V] = svd(A);
	catch e
		disp('Problem with intersection of planes');
		V = eye(3);
		S = eye(3);
	end
	L = plucker_dm_to_matrix(V(:,3),zeros(3,1));
	err = S(3,3);
	G = zeros(size(M));
	%figure;
	%plot3([G(1,:);M(1,:)],[G(2,:);M(2,:)],[G(3,:);M(3,:)],'-b');
	%hold on;
	%plot3([0,V(1,1)],[0,V(2,3)],[0,V(3,3)],'-r');
	%axis equal;
end
