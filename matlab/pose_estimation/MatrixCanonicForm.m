function [InvW, Qcanonic] = MatrixCanonicForm( Q )
% Matrix InW such that
% Qcanonic = transpose(InvW)*Q*InvW is diagonal

par         = [ Q(1,1); 2*Q(1,2); Q(2,2); 2*Q(1,3); 2*Q(2,3); Q(3,3)];
thetarad    = 0.5*atan2(par(2),par(1) - par(3));
cost        = cos(thetarad);
sint        = sin(thetarad);
sin_squared = sint * sint;
cos_squared = cost * cost;
cos_sin     = sint * cost;

Ao          = par(6);
Au          = par(4) * cost + par(5) * sint;
Av          = -par(4) * sint + par(5) * cost;
Auu         = par(1) * cos_squared + par(3) * sin_squared + par(2) * cos_sin;
Avv         = par(1) * sin_squared + par(3) * cos_squared - par(2) * cos_sin;

tuCentre    = - Au/(2*Auu);
tvCentre    = - Av/(2*Avv);

uCentre     = tuCentre * cost - tvCentre * sint;
vCentre     = tuCentre * sint + tvCentre * cost;
   
ctrQ        = [uCentre;vCentre];;
T           = [ eye(2,2), -ctrQ; 0, 0, 1];
InvT        = [ eye(2,2), ctrQ; 0, 0, 1];
R           = [ cost,-sint,0; sint,cost,0; 0, 0, 1 ];
InvW        = InvT*R;
Qcanonic    = transpose(InvW)*Q*InvW;

% Check in case we are not well oriented
if abs(Qcanonic(2,2)) > abs(Qcanonic(1,1))
	R = [ -cost,sint,0; -sint,-cost,0; 0, 0, 1 ];
	InvW = InvT*R;
	Qcanonic = transpose(InvW)*Q*InvW;
end

end
