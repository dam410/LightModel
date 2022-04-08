% This function estimate the possible plane perspectivity from intrinsic matrix
%	and local affine correspondence
function [h_vects] = affine_to_perspectivity(K,H_aff,pt_barycenter)

	% Augment barycenter vector
	X_b = [pt_barycenter;1];

	% Calculate image of this point with affine transform
	X_im_aff = H_aff*X_b;
	r_c_aff = norm(X_im_aff(1:2));

	% X_c is the barycenter of the 2D points in image coordinates
	K_inv = inv(K);

	% Calculate the barycenter point in normalized image coordinates
	X_normalized = K_inv*X_b;

	% Calculate the equivalent affine transformation that map 
	%	[0,0,0] centered at a circle point (not the center) to barycenter 
	%	in normalized image coordinates
	J_H = K_inv*inv(H_aff)*[eye(2),H_aff(1:2,:)*X_b;[0,0,1]];

	% Apply the IPPE algorithm
	[R1,R2,gamma] = IPPE_dec(X_normalized(1:2),J_H(1:2,1:2));
	
	% Calculate both possible homographies
	t_1 = 1/gamma*X_normalized-[R1(:,1:2),zeros(3,1)]*H_aff*X_b;
	t_2 = 1/gamma*X_normalized-[R2(:,1:2),zeros(3,1)]*H_aff*X_b;
	H_1 = K*[R1(:,1:2),t_1];
	H_2 = K*[R2(:,1:2),t_2];


	% Check if the homography keep the same radius for the barycenter point
	X_im_1 = inv(H_1)*X_b;
	X_im_1 = X_im_1./X_im_1(3);
	r_c_1 = norm(X_im_1(1:2));


	% Calculate the inverse homographies and vectorize it
	inv_H_1 = inv(H_1);	
	inv_H_2 = inv(H_2);	
	inv_H_1 = inv_H_1./inv_H_1(3,3);
	inv_H_2 = inv_H_2./inv_H_2(3,3);
	h_vects = transpose([inv_H_1(1:8);inv_H_2(1:8)]);

	H_1
	
end
