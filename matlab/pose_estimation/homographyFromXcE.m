% This function which intend to replace TransformationFromP (hard to understand and to use)

% Function homographyFromXcE
% brief : This function computes an homography from the image of a circle, its center and its rotation angle.
% Input Parameters :
%	Xc : a 2x1 vector which is the image of the center of the circle in the image coordinates
%	E : a 3x3 symmetric matrix which is the image of the circle in the image coordinaes
%	angle : a rotation angle
%
% Output Parameters :
%	invH : Homography mapping a point X in image coordinate to the rectified coordinate where the circle
% 	is the unit circle.
%	d_invH_d_xc : Derivatives of the inverse of the matrix along the first component of the image center coordinates
%	d_invH_d_yc : Derivatives of the inverse of the matrix along the first component of the image center coordinates
%
% 
function [invH,d_invH_d_xc,d_invH_d_yc] = homographyFromXcE(Xc,E,angle)

[InvW_,Qcanonic_] = MatrixCanonicForm(E);
W = inv(InvW_);
% We put Qcanonic in the special form
%[ 1/a^2,     0,  0]
%[     0, 1/b^2,  0]
%[     0,     0, -1]
Qcanonic_ = -Qcanonic_/Qcanonic_(3,3);
l3 = Qcanonic_(3,3);
l2 = Qcanonic_(2,2);
l1 = Qcanonic_(1,1);


Xc_ = W*[Xc;1];
Xc_ = Xc_/Xc_(3);
xc = Xc_(1);
yc = Xc_(2);

HcenterA = ...
[l3, l2*xc*yc, -l3*xc;...
0, - l1*xc^2 - l3, -l3*yc;...
-l1*xc, l2*yc, -l3];
HcenterB = diag([ ((l2*l3/l1*(l1*xc^2 + l2*yc^2 + l3)))^(1/2);l3;(-l2*(l1*xc^2 + l3))^(1/2)]);
Hcenter = HcenterA*HcenterB;
invHcenter =  inv(Hcenter);
R = [cos(angle),-sin(angle),0;sin(angle),cos(angle),0;0,0,1];
H = InvW_*Hcenter*R;
invH = inv(H);
% Test if point is outside the ellipse : in this case H will be complex
if ([Xc;1]'*E*[Xc;1]) > 0
	invH = eye(3);
end

if (nargout > 1)
	invHcenterA = inv(HcenterA);
	invHcenterB = inv(HcenterB);
	d_invHa_d_xc = ...
	[-(2*l1*xc)/(l1*xc^2 + l3)^2,...
	0,...
	(2*l1*xc^2)/(l1*xc^2 + l3)^2 - 1/(l1*xc^2 + l3);...
	(l1*yc)/(l1^2*xc^4 + 2*l1*l3*xc^2 + l2*l1*xc^2*yc^2 + l3^2 + l2*l3*yc^2) ...
	- (l1*xc*yc*(4*l1^2*xc^3 + 2*l2*l1*xc*yc^2 + 4*l3*l1*xc))/(l1^2*xc^4 + 2*l1*l3*xc^2 ...
	+ l2*l1*xc^2*yc^2 + l3^2 + l2*l3*yc^2)^2,...
	(2*l1*xc)/(l1*xc^2 + l2*yc^2 + l3)^2,...
	-(l3*yc*(4*l1^2*xc^3 + 2*l2*l1*xc*yc^2 + 4*l3*l1*xc))/(l1^2*xc^4 + 2*l1*l3*xc^2 ...
	+ l2*l1*xc^2*yc^2 + l3^2 + l2*l3*yc^2)^2;...
	(2*l1^2*l3*xc^2)/(l3^2 + l1*l3*xc^2 + l2*l3*yc^2)^2 - l1/(l3^2 + l1*l3*xc^2 + l2*l3*yc^2),...
	(2*l1*l2*l3*xc*yc)/(l3^2 + l1*l3*xc^2 + l2*l3*yc^2)^2,...
	(2*l1*xc)/(l1*xc^2 + l2*yc^2 + l3)^2];

	d_invHa_d_yc = ...
	[0,...
	0,...
	0;...
	(l1*xc)/(l1^2*xc^4 + 2*l1*l3*xc^2 + l2*l1*xc^2*yc^2 + l3^2 + l2*l3*yc^2) ...
	- (l1*xc*yc*(2*l1*l2*yc*xc^2 + 2*l2*l3*yc))/(l1^2*xc^4 + 2*l1*l3*xc^2 ...
	+ l2*l1*xc^2*yc^2 + l3^2 + l2*l3*yc^2)^2,...
	(2*l2*yc)/(l1*xc^2 + l2*yc^2 + l3)^2,...
	l3/(l1^2*xc^4 + 2*l1*l3*xc^2 + l2*l1*xc^2*yc^2 + l3^2 + l2*l3*yc^2) ...
	- (l3*yc*(2*l1*l2*yc*xc^2 + 2*l2*l3*yc))/(l1^2*xc^4 + 2*l1*l3*xc^2 ...
	+ l2*l1*xc^2*yc^2 + l3^2 + l2*l3*yc^2)^2;...
	(2*l1*l2*l3*xc*yc)/(l3^2 + l1*l3*xc^2 + l2*l3*yc^2)^2,...
	(2*l2^2*l3*yc^2)/(l3^2 + l1*l3*xc^2 + l2*l3*yc^2)^2 - l2/(l3^2 + l1*l3*xc^2 + l2*l3*yc^2),...
	(2*l2*yc)/(l1*xc^2 + l2*yc^2 + l3)^2];

	d_invHb_d_xc = ...
	[ -(l2*l3*xc)/((l2*l3*(l1*xc^2 + l2*yc^2 + l3))/l1)^(3/2), 0,0;...
	0, 0,0;...
	0, 0, (l1*l2*xc)/(-l2*(l1*xc^2 + l3))^(3/2)];

	d_invHb_d_yc = ...
	[ -(l2^2*l3*yc)/(l1*((l2*l3*(l1*xc^2 + l2*yc^2 + l3))/l1)^(3/2)), 0, 0;...
	0, 0, 0;...
	0, 0, 0];

	d_invH_d_xc = R'*(d_invHb_d_xc*invHcenterA + invHcenterB*d_invHa_d_xc)*W;
	d_invH_d_yc = R'*(d_invHb_d_yc*invHcenterA + invHcenterB*d_invHa_d_yc)*W;

end