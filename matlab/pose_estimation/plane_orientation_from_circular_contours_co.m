% Function to estimate the orientation of the scene plane using one isocontour
%	We distinguih K (only focal, ratio and sensor size) with K_im (take into account also the resolution in pixel
%		and the axial inversion)
% Specific for the colocalized PLS and optical center
function [X,N] = plane_orientation_from_circular_contours_co(K_im,T_cam,circular_isocontours,pt_vis);
	% Calculate orientation of the plane from circular isocontours
	nb_iso = 1;%length(circular_isocontours);
	N_vectors = zeros(3,nb_iso);
	X_vectors = zeros(3,nb_iso);
	% Check if a visible point is provided for convention on normal orientation
        if nargin < 4
                pt_vis = [];
        end
	for i=1:nb_iso
		ell = transpose_ellipse(circular_isocontours{i});
		E = param2ellipse(ell);
		% Direct method adapt for colocalized source/camera
		Q = K_im'*E*K_im;
		N_ = vanishing_line_paper_co(Q);
		N = N_(:,1);
		X = K_im*inv(Q)*N;
		X = X/X(3);
		N_vectors(:,i) = N;
		X_vectors(:,i) = X;
	end
	N = mean(N_vectors,2);
	N = N/norm(N);
	X = mean(X_vectors,2);
	% Apply the convention for the normal
	N = convention_normals(N,X,K_im,pt_vis);
end
