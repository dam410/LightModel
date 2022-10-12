% Special case 7 where the PLS is colocalised with camera
%	For this case as we cannot solve the ambiguity by using the intersection because
%	multiple planes does not share common variables.
%	We might as well use directly IPPE from a visible point
function [dssp,psp,ps] = case_7_co_ippe_ps(data)
	% Known parameters
	ps = {[0;0;0]};
	K = data.K;
	% Estimation (Using the homography directly and one special visible point)
	n_p = length(data.homography);

	dssp = cell(1,n_p);
	psp = cell(1,n_p);

	for i_p = 1:n_p
		pt_vis = get_visible_point_from_data(data,i_p);
		% The visible point coordinate
		w_normalised = inv(K)*[pt_vis(1);pt_vis(2);1]
		% Decompose the homography using IPPE but on the visible point
		%	Only extract the homography
		H_p = data.homography{i_p}{1};
		H_normalised = H_p*K;
		% Image of the visible point on the rectified plane (isophotes are circles)
		t = H_normalised*w_normalised;
		t = t/t(3);
		% Homography to change rectified plane coordinates so that 
		%	visible point moves to [0;0;1] once rectified
		H_center = [1,0,-t(1);0,1,-t(2);0,0,1]*H_normalised;
		inv_H_center = inv(H_center);
		H = inv_H_center./inv_H_center(3,3);
		H*[0;0;1]
		%	Replace the visible point
		J = zeros(2,2);
		J(1,1) = H(1,1)-H(3,1)*H(1,3);
		J(1,2) = H(1,2)-H(3,2)*H(1,3);
		J(2,1) = H(2,1)-H(3,1)*H(2,3);
		J(2,2) = H(2,2)-H(3,2)*H(2,3);
		v = [H(1,3);H(2,3)];
		[R1,R2] = IPPE_dec(v,J);
		[R1,R2,gamma] = IPPE_dec(v,J);
		psp{i_p} = {real([R1(:,3);0.0]),real([R2(:,3);0.0])};
		dssp{i_p} = {Inf,Inf};
	end
end
