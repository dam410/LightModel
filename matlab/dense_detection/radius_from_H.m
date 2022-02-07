%function [r_pt,J_r_pt] = radius_from_H(pt_in,H)
function [r_pt,J_r_pt] = radius_from_H(pt_in,H)

	nb_pt = size(pt_in,1);
	% Project the points
	proj_pt = double([pt_in,ones(nb_pt,1)]*transpose(H));
	Qx = proj_pt(:,1);
	Qy = proj_pt(:,2);
	Qz = proj_pt(:,3);
	Qx2 = Qx.^2;
	Qy2 = Qy.^2;
	Qz2 = Qz.^2;
	% Calculate radius
	%Qx_Qz = (Qx./Qz).^2;
	%Qy_Qz = (Qy./Qz).^2;
	%r_pt = Qx_Qz+Qy_Qz;
	Qx2Qy2 = Qx2+Qy2;
	r_pt = Qx2Qy2./Qz2;
	%disp('Min and max in r_pt');
	%min(r_pt)
	%max(r_pt)

	% Calculate the Jacobian
	J_Qx = [pt_in(:,1),zeros(nb_pt,2),pt_in(:,2),zeros(nb_pt,2),ones(nb_pt,1),zeros(nb_pt,2)];
	J_Qy = [zeros(nb_pt,1),pt_in(:,1),zeros(nb_pt,2),pt_in(:,2),zeros(nb_pt,2),ones(nb_pt,1),zeros(nb_pt,1)];
	J_Qz = [zeros(nb_pt,2),pt_in(:,1),zeros(nb_pt,2),pt_in(:,2),zeros(nb_pt,2),ones(nb_pt,1)];
	%J_r_pt = 2*(Qx.*J_Qx+Qy.*Qy)$
	%J_r_pt = 2*Qx_Qz.*(Qz.*J_Qx-Qx.*J_Qz)./Qz.^2+2*Qy_Qz.*(Qz.*J_Qy-Qy.*J_Qz)./Qz.^2;
	J_Qx2 = 2*Qx.*J_Qx;
	%disp('Min and max in J_Qx');
	%min(J_Qx(:))
	%max(J_Qx(:))
	J_Qy2 = 2*Qy.*J_Qy;
	J_Qz2 = 2*Qz.*J_Qz;
	Qz2_inv = 1./Qz2;
	J_Qz2_inv = -J_Qz2./Qz2;
	%disp('Min and max in J_Qz2_inv');
	%min(J_Qz2_inv(:))
	%max(J_Qz2_inv(:))
	%J_r_pt = Qz2_inv.*J_Qx2+Qx.^2.*J_Qz2_inv+Qz2_inv.*J_Qy2+Qy.^2.*J_Qz2_inv;
	J_r_pt = Qz2_inv.*(J_Qx2+J_Qy2)+Qx2Qy2.*J_Qz2_inv;
end
