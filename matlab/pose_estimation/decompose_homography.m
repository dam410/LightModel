function [R,t,scale] = decompose_homography(invK,invH)
	H = invK*inv(invH);
	[U,S,V] = svd(H(:,1:2));
	R_12 = U(:,1:2)*transpose(V);
	N = cross(R_12(:,1),R_12(:,2));
	R = [R_12,N];
	lambda = trace(transpose(R_12)*H(:,1:2))/trace(transpose(H(:,1:2))*H(:,1:2));
	t = H*[0;0;lambda];
end
