% This is a special implementation of monotonic cubic hermite spline
function [pp,J_pp_mat,J_b] = monotone_spline(x,y)
	n = length(x);
	% Calculate the tangents (m) at each points (x_i,y_i)

	%	Calculate the secant lines
	D = x(2:n)-x(1:(n-1));
	A = y(2:n)-y(1:(n-1));
	delta = A./D;
	%	Calculate the tangent 
	m_pre = [delta(1);(delta(1:(n-2))+delta(2:(n-1)))/2;delta(n-1)];
	%%%	Calculate (alpha,beta) vector used to check the strict monotony of 
	%%%	the interpolant function
	if nargin > 1
		J_x = [eye(n),zeros(n)];
		J_y = [zeros(n),eye(n)];
		J_D = J_x(2:n,:)-J_x(1:(n-1),:);
		J_inv_D = -J_D./D.^2;
		J_A = J_y(2:n,:)-J_y(1:(n-1),:);
		J_delta = J_A./D+J_inv_D.*A;
		J_m_pre = [J_delta(1,:);(J_delta(1:(n-2),:)+J_delta(2:(n-1),:))/2;J_delta(n-1,:)];
		J_alpha = zeros(n-1,2*n);
		J_beta = zeros(n-1,2*n);
		J_pp_mat = zeros(n-1,4,2*n);
		J_b = J_x;
	end
	alpha = zeros(1,n-1);
	beta = zeros(1,n-1);
	m = m_pre;
	J_m = J_m_pre;
	tau_vect = zeros(1,n-1);
	J_tau_vect = zeros(n-1,2*n);
	for i=1:(n-1)
		alpha(i) = m(i)/delta(i);
		beta(i) = m(i+1)/delta(i);
		J_alpha(i,:) = 1/(delta(i)^2).*(delta(i)*J_m(i,:)-m(i)*J_delta(i,:));
		J_beta(i,:) = 1/(delta(i)^2).*(delta(i)*J_m(i+1,:)-m(i+1)*J_delta(i,:));
		tau = 3/sqrt(alpha(i)^2+beta(i)^2);
		J_tau = [-(3*alpha(i))/(alpha(i)^2 + beta(i)^2)^(3/2),-(3*beta(i))/(alpha(i)^2 + beta(i)^2)^(3/2)]...
			*[J_alpha(i,:);J_beta(i,:)];
		tau_vect(i) = tau;
		J_tau_vect(i,:) = J_tau;
		if tau < 1
			m(i) = tau*alpha(i)*delta(i);
			m(i+1) = tau*beta(i)*delta(i);
			if nargout>1
				J_m(i,:) = alpha(i)*delta(i)*J_tau+tau*delta(i)*J_alpha(i,:)+tau*alpha(i)*J_delta(i,:);
				J_m(i+1,:) = beta(i)*delta(i)*J_tau+tau*delta(i)*J_beta(i,:)+tau*beta(i)*J_delta(i,:);
			end
		end
	end

	% Hermite Cubic basis matrix
	H = [2,-3,0,1;1,-2,1,0;-2,3,0,0;1,-1,0,0];
	% Initialize pp matrix
	pp_mat = zeros(n-1,4);
	for i=1:(n-1)
		if nargout > 1
			[H_t,J_H_t] = apply_affine_to_3_poly(H,1/D(i),0,J_inv_D(i,:));
			J_pp_mat(i,:,:) = ...
			transpose(H_t(1,:))*J_y(i,:) + permute(y(i)*J_H_t(1,:,:),[2,3,1]) + ...
			m(i)*transpose(H_t(2,:))*J_D(i,:) + D(i)*transpose(H_t(2,:))*J_m(i,:) + permute(D(i)*m(i)*J_H_t(2,:,:),[2,3,1]) + ...
			transpose(H_t(3,:))*J_y(i+1,:) + permute(y(i+1)*J_H_t(3,:,:),[2,3,1]) + ...
			m(i+1)*transpose(H_t(4,:))*J_D(i,:) + D(i)*transpose(H_t(4,:))*J_m(i+1,:) + permute(D(i)*m(i+1)*J_H_t(4,:,:),[2,3,1]);


		else
			[H_t] = apply_affine_to_3_poly(H,1/D(i),0);
		end
		pp_mat(i,:) = y(i)*H_t(1,:)+D(i)*m(i)*H_t(2,:)+y(i+1)*H_t(3,:)+D(i)*m(i+1)*H_t(4,:);
	end
	% Create the matlab pp struct
	pp = struct();
	pp.form = 'pp';
	pp.breaks = transpose(x(:));
	pp.coefs = pp_mat;
	pp.pieces = n-1;
	pp.order = 4;
	pp.dim = 1;
end

% Calculate a new 3 degree polynomial after an affine application Q(X) = P(a.X+b)
% Note : P and Q can be matrix
%	For the Jacobian we consider P, constant and b=0
function [Q,J_Q] = apply_affine_to_3_poly(P,a,b,J_a)
	M = [   a^3,3*a^2*b,3*a*b^2,b^3;...
		0,a^2,2*a*b,b^2;...
		0,0,a,b;...
		0,0,0,1];
	Q = P*M;
	if nargin>3
		nb_var = size(J_a,2);
		diag_J_M = [3*a^2*J_a;2*a*J_a;J_a;zeros(1,nb_var)];
		J_Q = zeros(size(P,1),size(P,2),nb_var);
		for i=1:size(P,1)
			J_Q(i,:,:) = transpose(P(i,:)).*diag_J_M;
		end
	end
end

