% This function estimate the possible plane perspectivity from intrinsic matrix
%	and local affine correspondence
function [h_vects,pt_centers] = homography_to_perspectivity(K,H_p,pt_barycenter)


	% Calculate the augmented vector
	w = [pt_barycenter;1];
	%w = inv(H_p)*[0;0;1];

	% Find a specific homography so that the distance between
	%	circles center and pt_barycenter on the scene plane is 1.
	% Je ne comprends vraiment pas pourquoi, ca ne fonctionne pas en considérant la
	%	translation et pourquoi une rotation du plan dégrade telle le calcul de la rotation
	% Alors qu'il devrait s'agir de la même pose du plan
	%w_p = H_p*w;
	%w_p = w_p/w_p(3);
	%scaling = 1;%norm(w_p(1:2));
	%alpha = atan2(w_p(2),w_p(1));
	%H_shift = diag([1/scaling,1/scaling,1])*[cos(alpha),-sin(alpha),0;sin(alpha),cos(alpha),0;0,0,1]*...
	%	[eye(2),-w_p(1:2);0,0,1];
	%H_p = H_shift*H_p;
	
	%w = inv(H_p)*w_p;
	%w = w/w(3);

	% Calculate the translation vector
	t = -H_p*w;
	t = t/t(3);
	
	% Check the projection of the point at the barycenter
	w_proj = H_p*w;
	r_barycenter = (w_proj(1)^2+w_proj(2)^2)/w_proj(3)^2;

	% Calculate the normalized cooridnate vectors
	w_normalized = inv(K)*w;
	w_normalized = w_normalized/w_normalized(3);
	
	% Calculate the composed homography
	H = inv(K)*inv(H_p)*[eye(2),+t(1:2);0,0,1];

	% Check that homography is correct
	X_test = H*[0;0;1];
	X_test = X_test/X_test(3);
	%w_normalized(1:2)-X_test(1:2)

	% Calculate the jacobian
	[X_0,J_H] = homography_jacobian(H);
	X_0-w_normalized;
	%X_0
	%J_H

	% Calculate the rotation matrix with IPPE
	[R1,R2,gamma] = IPPE_dec(X_0(1:2),J_H(1:2,1:2));
	R1
	R2

	 % Calculate both possible homographies
        t_1 = 1/gamma*w_normalized;%-[R1(:,1:2),zeros(3,1)]*H_p*w;
        t_2 = 1/gamma*w_normalized;%-[R2(:,1:2),zeros(3,1)]*H_p*w;
        H_1 = K*[R1(:,1:2),t_1]*[eye(2),-t(1:2);0,0,1];
        H_2 = K*[R2(:,1:2),t_2]*[eye(2),-t(1:2);0,0,1];

	% Calculate the inverse homographies and vectorize it
	inv_H_1 = inv(H_1);
	inv_H_2 = inv(H_2);
	inv_H_1 = inv_H_1./inv_H_1(3,3);
	inv_H_2 = inv_H_2./inv_H_2(3,3);

	pt_center_1 = H_1*[0;0;1];
	pt_center_2 = H_2*[0;0;1];
	pt_centers = {pt_center_1/pt_center_1(3),pt_center_2/pt_center_2(3)};

	% Calculate the same image of w
	w_proj_1 = inv_H_1*w;
	w_proj_2 = inv_H_2*w;
	r_barycenter_1 = (w_proj_1(1)^2+w_proj_1(2)^2)/w_proj_1(3)^2;
	r_barycenter_2 = (w_proj_2(1)^2+w_proj_2(2)^2)/w_proj_2(3)^2;
	h_vects = transpose([inv_H_1(1:8);inv_H_2(1:8)]);

end
