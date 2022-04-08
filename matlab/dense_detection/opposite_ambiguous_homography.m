function [h_vect_1,h_vect_2] = opposite_ambiguous_homography(K,H_opt,r_opt)
	invH_opt = inv(H_opt);
	Xc3 = invH_opt*[0;0;1];
	Xc3 = Xc3/Xc3(3);
	r_avg = mean(sqrt(r_opt(r_opt>0)));
	E = transpose(H_opt)*diag([1,1,-r_avg^2])*H_opt;
	ell_avg = ellipse2param(E);
	[Xc1,Xc2,JXc1,JXc2,N,JN] = image_center(K,ell_avg);
	Xc1
	Xc2
	Xc3
	pause;
	[invH_1] = homographyFromXcE(Xc1,E,0);
	invH_1 = diag([r_avg,r_avg,1])*invH_1;
	invH_1 = invH_1./invH_1(3,3);
	[invH_2] = homographyFromXcE(Xc2,E,0);
	invH_2 = diag([r_avg,r_avg,1])*invH_2;
	invH_2 = invH_2./invH_2(3,3);
	h_vect_1 = transpose(invH_1(1:8));
	h_vect_2 = transpose(invH_2(1:8));
	invH_1*invH_opt
end

