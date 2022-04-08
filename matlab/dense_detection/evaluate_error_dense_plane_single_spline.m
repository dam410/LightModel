% Evaluate the error and its jacobian for one plane
function [err,J_err] = evaluate_error_dense_plane_single_spline(inv_K,pts,I_pt,I_vect,NB_R,x_param,display)
	nb_pt = length(pts);
	EPS = 1;
	% Calculate the Jacobian for Xs (Projection of the source on the plane)
%	x_param = [Sx,Sy,Sz,alpha,beta,d,r_vect]
	S = x_param(1:3);
	[N,J_N_partial] = angle_to_normal(x_param(4),x_param(5));
	d = x_param(6);
	% Calculate the back projection of the points on the plane
	P = inv_K*transpose([pts,ones(nb_pt,1)]);
	% Calculate intersection with the plane
	interm = transpose(P)*N;
	inv_dist = 1./interm;
	%J_inv_dist_full = [zeros(nb_pt,3),J_inv_dist,zeros(nb_pt,1+NB_R)];
	lambdas = -d*inv_dist;
	P_scaled = transpose(lambdas).*P;
	% Calculate the value of the new parameter (no longer r distance on plane) 
	%	but u = cos(N,(P-S)/||P-S||)/(||P-S||^2 + EPS)
	%	which take into account the distance to the source (Global optimization possible)
	r_pt = sqrt((P_scaled(1,:)-S(1)).^2 + (P_scaled(2,:)-S(2)).^2 + (P_scaled(3,:)-S(3)).^2);
	cos_thetas = dot((P_scaled-S)./r_pt,N*ones(1,nb_pt));
	u_pt = transpose(cos_thetas./(r_pt.^2+ EPS));
	%figure, plot(r_pt);
	% Calculate the spline 
	% 	Import the r_vect vectors which is now different from directly the values in x_param
	J_r_vect = tril(ones(NB_R));
	a_vect = x_param(7:(NB_R+6));
	r_vect = J_r_vect*a_vect;
	if nargout > 1
		% Calculate the Jacobian
		J_N = [zeros(3),J_N_partial,zeros(3,1)];
		J_S = [eye(3),zeros(3,3)];
		J_d = zeros(1,6);
		J_d(6) = 1;
		J_interm = transpose(P)*J_N;
		J_inv_dist = (-1./interm.^2).*J_interm;
		J_lambdas = -d*J_inv_dist-inv_dist*J_d;
		% Calculate the distance between P and S
		%	Y_i = lambda_i.P(:,i) - S	
		%	J_yi = P(:,i).*J_lambda_i-J_S
		% 	r_i = ||Y_i||^2 = transpose(Y_i)*Y_i
		%	J_r_i = 2 transpose(Y_i)*J_yi
		%	      =   2*lambda_i transpose(P(:,i)) P(:,i) J_lambda_i (1)
		%		- 2*lambda_i transpose(P(:,i)) J_S              (2)
		%		- 2*transpose(S)*(P(:,i).*J_lambda_i)           (3)
		%		+ 2*transpose(S)*J_S                           (4)
		% All the term can be calculated separately
		J_r_pt_2 = 2*lambdas.*transpose(sum(P.*P,1)).*J_lambdas;% (1)
		J_r_pt_2 = J_r_pt_2 - 2*lambdas.*transpose(P)*J_S;% (2)
		J_r_pt_2 = J_r_pt_2 + 2*transpose(S)*J_S;% (4)
		J_r_pt_2 = J_r_pt_2 - 2*(transpose(P)*S).*J_lambdas;% (3)
		J_r_pt = 1./transpose(2.*r_pt).*J_r_pt_2;
		% Calculate the cosinus between normal and light source direction
		%	Y_i = lambda_i.P(:,i) - S
		%	J_yi =  P(:,i).*J_lambda_i-J_S
		%	V_i = Y_i/r_pt
		%	J_vi = J_yi/r_pt-Y_i.*J_r_pt/r_pt;
		%	cos_i = transpose(N)*V_i
		%	J_cosi = transpose(J_vi)*N + transpose(V_i)*J_N
		J_cost = transpose((P_scaled-S)./r_pt)*J_N;
		J_cost = J_cost + ((transpose(P)*N)./transpose(r_pt)).*J_lambdas;
		J_cost = J_cost -(transpose(N)*J_S)./transpose(r_pt);
		J_cost = J_cost-transpose((transpose(N)*(P_scaled-S))./r_pt.^2).*J_r_pt;
		J_u_pt = 1./transpose(r_pt.^2+EPS).*J_cost - transpose(cos_thetas./(r_pt.^2+EPS).^2).*J_r_pt_2;
		%J_u_pt = J_r_pt;
	end
	%u_pt = transpose(r_pt);
	[pp,J_pp_mat,J_breaks] = monotone_spline(r_vect,transpose(I_vect));
	if nargout > 1
		J_pp_mat_x = multiplication_tensor_matrix(J_pp_mat(:,:,1:NB_R),J_r_vect);
		J_breaks_x = J_breaks(:,1:NB_R)*J_r_vect;
		[I_pt_proj,J_I_pt_proj] = ppval_derivative_independent(pp,J_pp_mat_x,J_breaks_x,u_pt,J_u_pt);
		J_err = J_I_pt_proj(:,[(NB_R+1):(NB_R+6),1:NB_R]);
		err = I_pt_proj-I_pt;
	else
		I_pt_proj = ppval(pp,u_pt);
		err = I_pt_proj-I_pt;
	end
	if display==1
		display_spline_function(r_vect,u_pt,I_vect,I_pt,'Curve radius/intensity');
	end
end
