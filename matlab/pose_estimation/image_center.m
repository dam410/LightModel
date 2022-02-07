function [Xc1,Xc2,JXc1,JXc2,N,JN] = image_center(K,param_ellipse)

% Conversion of the parameter space of ellipse and jacobian of E
E = param2ellipse(param_ellipse);
inv_K = inv(K);
invE = inv(E);
% Put ellipse in a special form using intrinsic parameters
Q = K'*E*K;

% Compute the vanishing line with the new code
if nargout > 2 
	[L1,L2,J_L1,J_L2] = search_vanishing_line(Q);
else
	[L1,L2] = search_vanishing_line_no_jacobian(Q);
end
L = [L1,L2,zeros(3,1)];
X = K*inv(Q)*L;
% Backproject in 2D plane point
[Xc1,grad_X1_plane] = backproj(X(:,1));
[Xc2,grad_X2_plane] = backproj(X(:,2));

global inv_jacobian_sym;


% Compute the Jacobian if necessary
if (nargout > 2)
	[E,Je] = convertE(param_ellipse);
	JE = [zeros(9,4),Je];
	% Compute Jacobian of Q = K^(t)*E*K 
	JQ = zeros(9,9);
	% With the tools on this file, we only have to do multiplication
	JK = zeros(9,9);
	JK(1,1) = 1;
	JK(5,2) = 1;
	JK(7,3) = 1;
	JK(8,4) = 1;
	[K_t,JK_t] = matrix_trans(K,JK);
	[EK,JEK] = mult_mat(E,JE,K,JK);
	[Q,JQ] = mult_mat(K_t,JK_t,EK,JEK);

	% Compute Jacobian of E^(-1) use symbolic expression evaluate at startup
	Jinv_E = inv_jacobian_sym(E(1,1),E(1,2),E(1,3),E(2,1),E(2,2),E(2,3),E(3,1),E(3,2),E(3,3));
	Jinv_Q = inv_jacobian_sym(Q(1,1),Q(1,2),Q(1,3),Q(2,1),Q(2,2),Q(2,3),Q(3,1),Q(3,2),Q(3,3));
	invQ = inv(Q);
	JinvE = [zeros(9,4),Jinv_E*Je];
	JinvQ = Jinv_Q*JQ;
	JK = zeros(9,9);
	JK(1,1) = 1;		
	JK(5,2) = 1;
	JK(7,3) = 1;
	JK(8,4) = 1;
	% Compute Jacobian of K^(-t)
	[invK_t,JinvK_t] = inv_K_t(K);	
	% Compute jacobian of L in input parameters format
	J_L = [J_L1;J_L2;zeros(3,9)];
	JL = J_L*JQ;
	% Use matrix mutliplication with Jacobian
	[QL,JQL] = mult_mat(invQ,JinvQ,L,JL);
	[X,JX] = mult_mat(K,JK,QL,JQL);
	% Use Jacobian of backprojection
	JXc1 = grad_X1_plane*JX(1:3,:);
	JXc2 = grad_X2_plane*JX(4:6,:);

	% Compute also the normals and their derivative
	if nargout > 4
		% We use formula n = K l_{\infty}
		% Normalisation of N may reduce the value of the derivative and bring
		%	moe stability for ellipsoid computation
		n_front = [0;0;1];
		[n1,Jn1] = normalisation(L(:,1),JL(1:3,:));
		if n1'*n_front < 0
			n1 = -n1;
			Jn1 = -Jn1;
		end
		[n2,Jn2] = normalisation(L(:,2),JL(4:6,:));
		if n2'*n_front < 0
			n2 = -n2;
			Jn2 = -Jn2;
		end
		N = [n1,n2,zeros(3,1)];
		JN = zeros(9,9);
		JN(1:3,:) = Jn1;
		JN(4:6,:) = Jn2;
	end
end

end

function [L1,L2] = search_vanishing_line_no_jacobian(Q)

	b        = Q(1,1) + Q(2,2) + Q(3,3);
	c        = Q(1,1)*Q(2,2) + Q(1,1)*Q(3,3) + Q(2,2)*Q(3,3)-Q(1,2)^2 -Q(1,3)^2-Q(2,3)^2;
	d        = Q(1,1)*Q(2,3)^2 + Q(2,2)*Q(1,3)^2 + Q(3,3)*Q(1,2)^2 -Q(1,1)*Q(2,2)*Q(3,3)-2*Q(1,2)*Q(1,3)*Q(2,3);

	p        = b^2-3*c; % should not be zero or negative
	q        = 2*b^3-9*b*c-27*d;
	D        = acos(q/2/sqrt(p^3))
	if abs(D-pi) < 1e-7
		D = pi;
	end

	lambda   = 1/3*( b + 2*sqrt(p)*cos((D-2*pi)/3) );
	I3x3     = eye(3);
	X        = [0,Q(2:end)-lambda*I3x3(2:end)];
	X(1)     = (X(9)*X(4)^2 - 2*X(8)*X(4)*X(7) + X(5)*X(7)^2);
	s        = (- X(8)^2 + X(5)*X(9));
	X(2:end) = s*X(2:end);
	X        = reshape(X,3,3);

	v        = [(- X(2,3)^2 + X(2,2)*X(3,3)); (X(1,3)*X(2,3) - X(1,2)*X(3,3));(X(1,2)*X(2,3) - X(1,3)*X(2,2))];
	v        = v / sqrt(sum(v.*v));
	r1       = [-(X(1,3)*X(2,3) - X(1,2)*X(3,3));(- X(2,3)^2 + X(2,2)*X(3,3));0];
	r1       = r1 / sqrt(sum(r1.*r1));
	r2       = [-(X(1,2)*X(2,3) - X(1,3)*X(2,2))*(- X(2,3)^2 + X(2,2)*X(3,3)); ...
			-(X(1,2)*X(2,3) - X(1,3)*X(2,2))*(X(1,3)*X(2,3) - X(1,2)*X(3,3)); ...
			(X(1,3)*X(2,3) - X(1,2)*X(3,3))^2 + (- X(2,3)^2 + X(2,2)*X(3,3))^2];
	r2       = r2 / sqrt(sum(r2.*r2));
	R        = [ r1,r2,v];

	C        = transpose(R)*X*R;
	v        = v / sqrt(sum(v.*v));

	mu1      = C(1,1)/2 + C(2,2)/2 - (C(1,1)^2 - 2*C(1,1)*C(2,2) + 4*C(1,2)^2 + C(2,2)^2)^(1/2)/2;
	mu2      = C(1,1)/2 + C(2,2)/2 + (C(1,1)^2 - 2*C(1,1)*C(2,2) + 4*C(1,2)^2 + C(2,2)^2)^(1/2)/2;
	w1       = [ (C(1,1)/2 + C(2,2)/2 - (C(1,1)^2 - 2*C(1,1)*C(2,2) + 4*C(1,2)^2 + C(2,2)^2)^(1/2)/2) - C(2,2); C(1,2)  ];
	w2       = [ (C(1,1)/2 + C(2,2)/2 + (C(1,1)^2 - 2*C(1,1)*C(2,2) + 4*C(1,2)^2 + C(2,2)^2)^(1/2)/2) - C(2,2); C(1,2)  ];
	w1       = w1 / sqrt(sum(w1.*w1));
	w2       = w2 / sqrt(sum(w2.*w2));

	L1 =R*[w1*sqrt(-mu1)+w2*sqrt(mu2);0];
	L2 =R*[w1*sqrt(-mu1)-w2*sqrt(mu2);0];
end

function [L1,L2,JL1,JL2] = search_vanishing_line(Q)

% Read the matrix variable exactly as reshape, column by column
b        = Q(1,1) + Q(2,2) + Q(3,3);
Jb = [1 0 0 0 1 0 0 0 1];

c        = Q(1,1)*Q(2,2) + Q(1,1)*Q(3,3) + Q(2,2)*Q(3,3)-Q(1,2)^2 -Q(1,3)^2-Q(2,3)^2;
Jc = [Q(2,2)+Q(3,3), 0, 0, -2*Q(1,2), Q(1,1)+Q(3,3), 0, -2*Q(1,3), -2*Q(2,3), Q(1,1)+Q(2,2)];

d        = Q(1,1)*Q(2,3)^2 + Q(2,2)*Q(1,3)^2 + Q(3,3)*Q(1,2)^2 ...
-Q(1,1)*Q(2,2)*Q(3,3)-2*Q(1,2)*Q(1,3)*Q(2,3);

Jd = [Q(2,3)^2-Q(2,2)*Q(3,3), 0, 0, 2*Q(3,3)*Q(1,2)-2*Q(1,3)*Q(2,3), Q(1,3)^2-Q(1,1)*Q(3,3), 0, 2*Q(2,2)*Q(1,3)-2*Q(1,2)*Q(2,3),2*Q(1,1)*Q(2,3)-2*Q(1,2)*Q(1,3), Q(1,2)^2-Q(1,1)*Q(2,2)];


p        = b^2-3*c; % should not be zero or negative
Jp = 2*b*Jb-3*Jc;

q        = 2*b^3-9*b*c-27*d;
Jq = 2*3*b^2*Jb-9*(b*Jc+c*Jb)-27*Jd;

% D value can become complex if the ellipse is already circular, e.g. circle axis passing trough camera optical center
%	taking real component fix the problem but make the jacobian undefined (sqrt(0))
D        = real(acos(q/2/sqrt(p^3)));
JD = real(-(p^(-3/2)*Jq-3/2*q*p^(-5/2)*Jp)/(2*sqrt(1-q^2*p^(-3)/4)));

lambda   = 1/3*( b + 2*sqrt(p)*cos((D-2*pi)/3) );
Jlambda = 1/3*( Jb + Jp/sqrt(p)*cos((D-2*pi)/3) - 2/3*sqrt(p)*sin((D-2*pi)/3)*JD );

I3x3     = eye(3);
X        = [0,Q(2:end)-lambda*I3x3(2:end)];
JX = eye(9,9) - I3x3(:)*Jlambda;

X(1)     = (X(9)*X(4)^2 - 2*X(8)*X(4)*X(7) + X(5)*X(7)^2);
JX(1,:) = X(4)^2*JX(9,:)+2*X(9)*X(4)*JX(4,:)-2*(X(8)*X(4)*JX(7,:)+X(8)*X(7)*JX(4,:)+X(4)*X(7)*JX(8,:))+X(7)^2*JX(5,:)+2*X(5)*X(7)*JX(7,:);

s        = (- X(8)^2 + X(5)*X(9));
Js = -2*X(8)*JX(8,:) + X(5)*JX(9,:) + X(9)*JX(5,:);

X_save(2:9) = X(2:9);
X(2:end) = s*X(2:end);
JX(2:9,:) = s*JX(2:9,:) + X_save(2:9)'*Js;

X        = reshape(X,3,3);

v        = [(- X(2,3)^2 + X(2,2)*X(3,3)); (X(1,3)*X(2,3) - X(1,2)*X(3,3));(X(1,2)*X(2,3) - X(1,3)*X(2,2))];
% The translation in line coordinates give
% v = [ -X8^2 + X5*X9 ; X7*X8 - X4*X9 ; X4*X8 - X7*X5 ]
Jv = [-2*X(8)*JX(8,:) + X(5)*JX(9,:) + X(9)*JX(5,:) ; X(7)*JX(8,:) + X(8)*JX(7,:) - X(4)*JX(9,:) - X(9)*JX(4,:) ;...
	X(4)*JX(8,:) + X(8)*JX(4,:) - X(7)*JX(5,:) - X(5)*JX(7,:) ];
%v        = v / sqrt(sum(v.*v));
[v,Jv] = normalisation(v,Jv);
r1       = [-(X(1,3)*X(2,3) - X(1,2)*X(3,3));(- X(2,3)^2 + X(2,2)*X(3,3));0];
% r1 = [ - X7*X8 - X4*X9 ; -X8^2 + X5*X9 ; 0]
Jr1 = [ -X(7)*JX(8,:) - X(8)*JX(7,:) + X(4)*JX(9,:) + X(9)*JX(4,:) ;...
	-2*X(8)*JX(8,:) + X(5)*JX(9,:) + X(9)*JX(5,:) ;...
	zeros(1,9)];
%r1       = r1 / sqrt(sum(r1.*r1));
[r1,Jr1] = normalisation(r1,Jr1);
%r2       = [-(X(1,2)*X(2,3) - X(1,3)*X(2,2))*(- X(2,3)^2 + X(2,2)*X(3,3));...
%	-(X(1,2)*X(2,3) - X(1,3)*X(2,2))*(X(1,3)*X(2,3) - X(1,2)*X(3,3));...
%	(X(1,3)*X(2,3) - X(1,2)*X(3,3))^2 + (- X(2,3)^2 + X(2,2)*X(3,3))^2];
%r2 = [-(X4*X8 - X7*X5)*(-X8^2+X5*X9); -(X4*X8 - X7*X5)*(X7*X8 - X4*X9); (X7*X8 - X4*X9)^2 + (-X8^2 + X5*X9)^2];
% Simplification :
%	Ux = (X4*X8 - X7*X5)
Ux = X(4)*X(8)-X(7)*X(5);
JUx = X(4)*JX(8,:) + X(8)*JX(4,:) - X(7)*JX(5,:) - X(5)*JX(7,:);
%	Vx = (-X8^2+X5*X9)
Vx = -X(8)^2 + X(5)*X(9);
JVx = -2*X(8)*JX(8,:) + X(5)*JX(9,:) + X(9)*JX(5,:);
%	Wx = (X7*X8 - X4*X9)
Wx = X(7)*X(8) - X(4)*X(9);
JWx = X(7)*JX(8,:) + X(8)*JX(7,:) - X(4)*JX(9,:) - X(9)*JX(4,:);
% r2 = [ -Ux*Vx ; -Ux*Wx ; Wx^2 + Vx^2 ]
r2 = [-Ux*Vx ; -Ux*Wx ; Wx^2 + Vx^2];
Jr2 = [-Ux*JVx - Vx*JUx; -Ux*JWx - Wx*JUx ; 2*Wx*JWx + 2*Vx*JVx];
%r2       = r2 / sqrt(sum(r2.*r2));
[r2,Jr2] = normalisation(r2,Jr2);
R        = [r1,r2,v];
JR = [Jr1;Jr2;Jv];

%C        = transpose(R)*X*R
[C_,JC_] = mult_mat(X,JX,R,JR);

[R_t,JR_t] = matrix_trans(R,JR);
[C,JC] = mult_mat(R_t,JR_t,C_,JC_);

Delta = (C(1,1)^2 - 2*C(1,1)*C(2,2) + 4*C(1,2)^2 + C(2,2)^2)^(1/2)/2;

JDelta = 2*C(1,1)*JC(1,:) - 2*C(1,1)*JC(5,:) - 2*C(2,2)*JC(1,:) + 8*C(1,2)*JC(4,:) + 2*C(2,2)*JC(5,:);
JDelta = 1/sqrt(C(1,1)^2 - 2*C(1,1)*C(2,2) + 4*C(1,2)^2 + C(2,2)^2)*JDelta/4;

% mu1      = C(1,1)/2 + C(2,2)/2 - (C(1,1)^2 - 2*C(1,1)*C(2,2) + 4*C(1,2)^2 + C(2,2)^2)^(1/2)/2
mu1 = C(1,1)/2 + C(2,2)/2 - Delta;
Jmu1 = JC(1,:)/2+JC(5,:)/2 - JDelta;

% mu2      = C(1,1)/2 + C(2,2)/2 + (C(1,1)^2 - 2*C(1,1)*C(2,2) + 4*C(1,2)^2 + C(2,2)^2)^(1/2)/2
mu2 = C(1,1)/2 + C(2,2)/2 + Delta;
Jmu2 = JC(1,:)/2+JC(5,:)/2 + JDelta;

%w1       = [ (C(1,1)/2 + C(2,2)/2 - (C(1,1)^2 - 2*C(1,1)*C(2,2)...
%+ 4*C(1,2)^2 + C(2,2)^2)^(1/2)/2) - C(2,2); C(1,2)  ];
w1 = [mu1-C(2,2);C(1,2)];
Jw1 = [Jmu1-JC(5,:);JC(4,:)];

%w2       = [ (C(1,1)/2 + C(2,2)/2 + (C(1,1)^2 - 2*C(1,1)*C(2,2)...
%+ 4*C(1,2)^2 + C(2,2)^2)^(1/2)/2) - C(2,2); C(1,2)  ];
w2 = [mu2-C(2,2);C(1,2)];
Jw2 = [Jmu2-JC(5,:);JC(4,:)];

%w1       = w1 / sqrt(sum(w1.*w1));
%w2       = w2 / sqrt(sum(w2.*w2));
[w2,Jw2] = normalisation(w2,Jw2);
[w1,Jw1] = normalisation(w1,Jw1);

vectors = zeros(3,3);
vectors(1:2,1) = w1*sqrt(-mu1)+w2*sqrt(mu2);
vectors(1:2,2) = w1*sqrt(-mu1)-w2*sqrt(mu2);
Jvectors = [Jw1*sqrt(-mu1) - w1/2/sqrt(-mu1)*Jmu1 + Jw2*sqrt(mu2) + w2/2/sqrt(mu2)*Jmu2;zeros(1,9);Jw1*sqrt(-mu1) - w1/2/sqrt(-mu1)*Jmu1 - Jw2*sqrt(mu2) - w2/2/sqrt(mu2)*Jmu2 ;zeros(1,9);zeros(3,9)];
[L,JL] = mult_mat(R,JR,vectors,Jvectors);

%L1 =R*[w1*sqrt(-mu1)+w2*sqrt(mu2);0];
%[Jw1*sqrt(-mu1) - w1/2/sqrt(-mu1)*Jmu1 + Jw2*sqrt(mu2) + w2/2/sqrt(mu2)*Jmu2 ;zeros(1,9)]
%[L1,JL1] = mult_mat(R,JR,[w1*sqrt(-mu1)+w2*sqrt(mu2);0],[Jw1*sqrt(-mu1) - w1/2/sqrt(-mu1)*Jmu1 + Jw2*sqrt(mu2) + w2/2/sqrt(mu2)*Jmu2 ;zeros(1,9)]);

%L2 =R*[w1*sqrt(-mu1)-w2*sqrt(mu2);0]
%[L2,JL2] = mult_mat(R,JR,[w1*sqrt(-mu1)-w2*sqrt(mu2);0],[Jw1*sqrt(-mu1) - w1/2/sqrt(-mu1)*Jmu1 - Jw2*sqrt(mu2) - w2/2/sqrt(mu2)*Jmu2 ;zeros(1,9)]);

L1 = L(:,1);
L2 = L(:,2);
JL1 = JL(1:3,:);
JL2 = JL(4:6,:);


end

% Function use to compare calculated derivative with finite diff
function [Jdiff] = diffJ(fun,Q)
eps = 0.00001;
X0 = fun(Q);
size_in = size(Q);
size_out = size(X0);
n_in = prod(size_in);
n_out = prod(size_out);
Jdiff = zeros(n_out,n_in);
for i = 1:n_in
	Q_ = Q;
	Q_(i) = Q_(i) + eps;
	X1 = fun(Q_);
	D = (X1-X0)/eps;
	Jdiff(:,i) = D(:)';
end
end

% Compute the more generic normalisation of a vector
function [X_normed,JX_normed] = normalisation(X,JX)
norm_X = sqrt(X'*X);
Jnorm_X = X'*JX/norm_X;
JX_normed = 1/norm_X*JX - X*(Jnorm_X/(norm_X*norm_X));
X_normed = X/norm_X;
end

% Compute the Jacobian of the matrix multiplication X = A*B (nxn)
function [X,JX] = mult_mat(A,JA,B,JB)
[n_dim,~] = size(A);
[~,n_out] = size(JA);
X = A*B;
JX = trans_mult_left(A)*JB + trans_mult_right(B)'*JA;
end

% This transformation transform a matrix A in a matrix such that :
%	X = A*B
%	T(A)*B(:) = X(:)
function [A_mut] = trans_mult_left(A)
	[n_dim,~] = size(A);
	A_mut = zeros(n_dim*n_dim,n_dim*n_dim);
	for i = 1:n_dim
		A_mut((i-1)*n_dim+(1:n_dim),(i-1)*n_dim+(1:n_dim)) = A;
	end

end

% This transformation transform a matrix A in a matrix such that :
%	X = A*B
%	A(:)'*T(B) = X(:)
function [A_mut] = trans_mult_right(A)
	[n_dim,~] = size(A);
	Mask = repmat(eye(n_dim),n_dim,n_dim);
	A_mut = zeros(n_dim*n_dim,n_dim*n_dim);
	for i=0:(n_dim*n_dim-1)
		A_mut(n_dim*mod(i,n_dim)+(1:n_dim),n_dim*floor(i/n_dim)+(1:n_dim)) = A(i+1);
	end
	A_mut = A_mut.*Mask;
end


% Compute the transformation from a matrix and its jacobian to its jacobian of its transpose
function [A_t,JA_t] = matrix_trans(A,JA)
[n_dim,~] = size(A);
ind = 0:(n_dim*n_dim-1);
T = zeros(n_dim*n_dim);
T(ind*n_dim^2+3*mod(ind,3)+floor(ind/3)+1) = 1;
JA_t = T*JA;
A_t = A';
end


function [invK_t,JinvK_t] = inv_K_t(K)
	% We use the special form of K
	invK_t = [1/K(1,1), 0, 0;...
	0, 1/K(2,2), 0;...
	-K(1,3)/K(1,1), -K(2,3)/K(2,2), 1];
	JinvK_t = zeros(9,9);
	JinvK_t(1,1) = -1/(K(1,1)*K(1,1));
	JinvK_t(3,1) = K(1,3)/(K(1,1)*K(1,1));
	JinvK_t(5,2) = -1/(K(2,2)*K(2,2));
	JinvK_t(6,2) = K(2,3)/(K(2,2)*K(2,2));
	JinvK_t(3,3) = -1/K(2,2);
	JinvK_t(6,4) = -1/K(2,2);
end

% Function of back projection on plane z=1
function [X_plane,grad_X_plane] = backproj(X)
	X_plane = X(1:2)/X(3);
	grad_X_plane = [1/X(3), 0, -X(1)/(X(3)*X(3));...
	0, 1/X(3), -X(2)/(X(3)*X(3))];
end


