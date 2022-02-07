syms a11 a12 a13 a21 a22 a23 a31 a32 a33;
A_sym = [a11, a12, a13; a21, a22, a23; a31, a32, a33];
invA_sym = inv(A_sym);
invA_sym_lin = invA_sym(:);
Jinv_sym = jacobian(invA_sym_lin,[a11,a21,a31,a12,a22,a32,a13,a23,a33]);
global inv_jacobian_sym;
inv_jacobian_sym = matlabFunction(Jinv_sym);
