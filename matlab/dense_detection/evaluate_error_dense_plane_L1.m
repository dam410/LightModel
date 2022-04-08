% Evaluate the error and its gradient
function [norm_L1,gradient_norm_L1] = evaluate_error_dense_plane_L1(inv_K,pts,I_pt,I_vect,NB_R,x_param,display)
	nb_pt = length(pts);
	% Calculate the Jacobian for Xs (Projection of the source on the plane)
%	x_param = [Sx,Sy,Sz,alpha,beta,d,r_vect]
	S = x_param(1:3);
	[N,J_N_partial] = angle_to_normal(x_param(4),x_param(5));
	J_N = [zeros(3),J_N_partial,zeros(3,1)];
	J_S = [eye(3),zeros(3,3)];
	J_d = zeros(1,6);
	J_d(6) = 1;
	d = x_param(6);
	h =  -d-dot(S,N);
	J_h = -J_d-transpose(S)*J_N-transpose(N)*J_S;
	Xs = S+h*N;
	J_Xs = J_S + N.*J_h + h*J_N;
	% Calculate the back projection of the points on the plane
	P = inv_K*transpose([pts,ones(nb_pt,1)]);
	% Calculate intersection with the plane
	interm = transpose(P)*N;
	J_interm = transpose(P)*J_N;
	inv_dist = 1./interm;
	J_inv_dist = (-1./interm.^2).*J_interm;
	%J_inv_dist_full = [zeros(nb_pt,3),J_inv_dist,zeros(nb_pt,1+NB_R)];
	lambdas = -d*inv_dist;
	J_lambdas = -d*J_inv_dist-inv_dist*J_d;
	P_scaled = transpose(lambdas).*P;
	Y_i = lambdas(1)*P(:,1) - Xs;
	J_yi = P(:,1).*J_lambdas(1,:)-J_Xs;
	r_i = sum(Y_i.^2);
	J_r_i  = 2*transpose(Y_i)*J_yi;
	% Calculate the distance between P and Xs
	%	Y_i = lambda_i.P(:,i) - Xs	
	%	J_yi = P(:,i).*J_lambda_i-J_Xs
	% 	r_i = ||Y_i||^2 = transpose(Y_i)*Y_i
	%	J_r_i = 2 transpose(Y_i)*J_yi
	%	      =   2*lambda_i transpose(P(:,i)) P(:,i) J_lambda_i (1)
	%		- 2*lambda_i transpose(P(:,i)) J_Xs              (2)
	%		- 2*transpose(Xs)*(P(:,i).*J_lambda_i)           (3)
	%		+ 2*transpose(Xs)*J_Xs                           (4)
	% All the term can be calculated separately
	J_r_pt = 2*lambdas.*transpose(sum(P.*P,1)).*J_lambdas;% (1)
	J_r_pt = J_r_pt - 2*lambdas.*transpose(P)*J_Xs;% (2)
	J_r_pt = J_r_pt + 2*transpose(Xs)*J_Xs;% (4)
	J_r_pt = J_r_pt - 2*(transpose(P)*Xs).*J_lambdas;% (3)
	r_pt = transpose((P_scaled(1,:)-Xs(1)).^2+(P_scaled(2,:)-Xs(2)).^2+(P_scaled(3,:)-Xs(3)).^2);
	% Calculate the spline 
	% 	Import the r_vect vectors which is now different from directly the values in x_param
	J_r_vect = tril(ones(NB_R));
	x_vect = max(0,x_param(7:(NB_R+6)));
	r_vect = J_r_vect*x_vect;
	[pp,J_pp_mat,J_breaks] = monotone_spline(r_vect,transpose(I_vect));
	J_pp_mat_x = multiplication_tensor_matrix(J_pp_mat(:,:,1:NB_R),J_r_vect);
	J_breaks_x = J_breaks(:,1:NB_R)*J_r_vect;
	[I_pt_proj,J_I_pt_proj] = ppval_derivative_independent(pp,J_pp_mat_x,J_breaks_x,r_pt,J_r_pt);
	J_err = J_I_pt_proj(:,[(NB_R+1):(NB_R+6),1:NB_R]);
	err = I_pt_proj-I_pt;
	[norm_L1,gradient_norm_L1] = error_norm_L1(err,J_err);
	if display==1
		display_spline_function(r_vect,r_pt,I_vect,I_pt,'Curve radius/intensity');
	end
end
