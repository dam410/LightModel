% Apply a transformation to an ellipse given in parameters so that the ellipse can be displayed
%	accordingly in a transposed image
function [ell_t] = transpose_ellipse(ell)
	E = param2ellipse(ell);
	trans_coord = [0,1,0;1,0,0;0,0,1];
	E_t = trans_coord*E*trans_coord;
	ell_t = ellipse2param(E_t);
end
