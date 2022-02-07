function X = ellipse2param(ellipse)
% ELLIPSE2PARAM Extract the ellipse param from its matrix representation.
%
% SYNTAX
%   params = point_to_ellipse_distance(ellipse)
%
% INPUTS
%   - ellipse (3x3 matrix):
%       [  a, b/2, d/2;
%        b/2,   c, e/2;
%        d/2, e/2,   f ]
%
% OUTPUTS
%   - The parameters of the ellipse (x0, y0, a, b, theta)

a = ellipse(1, 1);
b = 2 * ellipse(1, 2);
c = ellipse(2, 2);
d = 2 * ellipse(1, 3);
e = 2 * ellipse(2, 3);
f = ellipse(3, 3);

par = [a b c d e f];

thetarad = 0.5 * atan2(par(2), par(1) - par(3));
cost = cos(thetarad);
sint = sin(thetarad);
sin_squared = sint .* sint;
cos_squared = cost .* cost;
cos_sin = sint .* cost;

Ao = par(6);
Au = par(4) .* cost + par(5) .* sint;
Av = - par(4) .* sint + par(5) .* cost;
Auu = par(1) .* cos_squared + par(3) .* sin_squared + par(2) .* cos_sin;
Avv = par(1) .* sin_squared + par(3) .* cos_squared - par(2) .* cos_sin;

if Auu == 0 || Avv == 0
   X = zeros(1, 5);
else
    tuCentre = - Au ./ (2 .* Auu);
    tvCentre = - Av ./ (2 .* Avv);
    wCentre  = Ao - Auu .* tuCentre .* tuCentre - Avv .* tvCentre .* tvCentre;

    uCentre = tuCentre .* cost - tvCentre .* sint;
    vCentre = tuCentre .* sint + tvCentre .* cost;

    Ru = -wCentre ./ Auu;
    Rv = -wCentre ./ Avv;

    Ru = sqrt(abs(Ru)).*sign(Ru);
    Rv = sqrt(abs(Rv)).*sign(Rv);

    X = [uCentre, vCentre, Ru, Rv, thetarad];
    % Garanti de la contrainte (Ru>Rv) et thetarad inclus
    % Assure that the ellipse is aligned on its major axis
    if X(3) < X(4)
    angle = mod(pi/2+X(5),pi);
        if angle>(pi/2)
            angle = angle-pi;
        end
        a = X(3);
        X(3) = X(4);
        X(4) = a;
        X(5) = angle;
    end
end
