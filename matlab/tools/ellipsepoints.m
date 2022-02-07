function pts = ellipsepoints(params,n,arc,T)
% function [h, pts] = ellipsepoints(params,n,T)
%          params (required) = [uc, vc, Ru, Rv, theta]; (theta in radians)
%          n      (required) = number of points
%          T      (optional) = homogeneous 3-by-3 transformation


if nargin < 3,
   arc = pi*2;
end    
if nargin < 4,
   T = eye(3,3);
end    

if params(3)*params(4) < 0
	t = linspace(arc/n,arc-(arc/n),n);
	x = params(3)./ cos(t);
	y = params(4) * tan(t);
else
	t   = linspace(0,arc,n);
	x = params(3) * cos(t);
	y = params(4) * sin(t);
end
nx  = x*cos(params(5))-y*sin(params(5)) + params(1); 
ny  = x*sin(params(5))+y*cos(params(5)) + params(2);

pts = T * augment([nx;ny]);
