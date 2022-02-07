% This function estimate the possible plane perspectivity from intrinsic matrix
%	and local affine correspondence
function [h_vects] = homography_to_ambiguous_perspectivity(invK,H_p)
	

	% Calculate the composed homography
	H = invK*inv(H_p);
	K = inv(invK);

	w = H*[0;0;1];
	w_normalized = w/norm(w);

	% Calculate the jacobian
	H = H./H(3,3);
	J = zeros(2,2);
	J(1,1) = H(1,1)-H(3,1)*H(1,3);
	J(1,2) = H(1,2)-H(3,2)*H(1,3);
	J(2,1) = H(2,1)-H(3,1)*H(2,3);
	J(2,2) = H(2,2)-H(3,2)*H(2,3);
	v = [H(1,3);H(2,3)];
	[R1,R2] = IPPE_dec(v,J);


	% Calculate the rotation matrix with IPPE
	[R1,R2,gamma] = IPPE_dec(v,J);

	 % Calculate both possible homographies
        %t_1 = 1/gamma*w_normalized;
        %t_2 = 1/gamma*w_normalized;
        H_1 = K*[R1(:,1:2),w_normalized];
        H_2 = K*[R2(:,1:2),w_normalized];
	
	inv_H_1 = inv(H_1);
	inv_H_1 = inv_H_1/inv_H_1(3,3);

	inv_H_2 = inv(H_2);
	inv_H_2 = inv_H_2/inv_H_2(3,3);

	% Calculate the homography in form of vector
	h_vects = [homography_to_vector_perspective(inv_H_1),...
		homography_to_vector_perspective(inv_H_2)];

end
