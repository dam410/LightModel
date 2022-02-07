function E = param2ellipse(param)
% PARAM2ELLIPSE Compute the ellipse matrix from its parameters.
%
% SYNTAX
%   E = param2ellipse(param)
%
% INPUTS
%   - param (5-vector): Parameters of the ellipse (x0, y0, a, b, theta).
%
% OUTPUTS
%   - E (matrix 3x3): the matrix of the ellipse (not normalized in any way).

if length(param)==6,
    error('For ellipse Coefficients to matrix: use coeffs2ellipse');
end

xc = param(1);
yc = param(2);
ax = param(3);
bx = param(4);
rho = param(5);

T = [cos(rho), -sin(rho), xc;
    sin(rho), cos(rho), yc;
    0, 0, 1];

invT = T \ eye(3);
E = transpose(invT) * diag([1/ax^2 ; 1/bx^2 ; -1]) * invT;

% (potential normalization)
% C = C/norm(C,'fro');
