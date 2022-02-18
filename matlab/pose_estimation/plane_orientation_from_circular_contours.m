% Function to estimate the orientation of the scene plane using one or two circular isocontours
%	We distinguih K (only focal, ratio and sensor size) with K_im (take into account also the resolution in pixel
%		and the axial inversion)
function [Xc_best,N_best] = plane_orientation_from_circular_contours(K_im,T_cam,circular_isocontours,pt_vis);
	% Calculate orientation of the plane from circular isocontours
	nb_iso = length(circular_isocontours);
	Xc = zeros(2*nb_iso,2);
	N = zeros(3*nb_iso,2);
	% Check if a visible point is provided for convention on normal orientation
	if nargin < 4 
		pt_vis = [];
	end
	for i=1:nb_iso
		ell = transpose_ellipse(circular_isocontours{i});
                E = param2ellipse(ell);
                % Direct method adapt for colocalized source/camera
                Q = K_im'*E*K_im;
		%[Xc1,Xc2,~,~,N_i] = image_center(K_im,ell);
		N_i = vanishing_line_paper_co(Q);
		% Now calculate the image of the brightest point
		if size(N_i,2)>1
			Xc1 = K_im*inv(Q)*N_i(:,1);
			Xc2 = K_im*inv(Q)*N_i(:,2);
			Xc1 = Xc1(1:2)/Xc1(3);
			Xc2 = Xc2(1:2)/Xc2(3);
		else
			Xc1 = K_im*inv(Q)*N_i(:,1);
			Xc2 = K_im*inv(Q)*N_i(:,1);
			Xc1 = Xc1(1:2)/Xc1(3);
			Xc2 = Xc2(1:2)/Xc2(3);
			N_i = [N_i,N_i];
		end
		Xc(2*i-1:2*i,1) = Xc1;
		Xc(2*i-1:2*i,2) = Xc2;
		%normals = transpose(inv(K_im))*N_i(:,1:2);
		normals = N_i(:,1:2);
		%N(3*i-2:3*i,:) = transpose(inv(K_im))*[normals(:,1)/norm(normals(:,1)),normals(:,2)/norm(normals(:,2))];
		N(3*i-2:3*i,:) = [normals(:,1)/norm(normals(:,1)),normals(:,2)/norm(normals(:,2))];
	end
	% Seems complicated but it actually only choose which ambiguous solution
	%	is correct by comparing the image center estimated position
	%	and select the ambiguous index that minimize it.
	if nb_iso>1
		comp_pos = @(x,y) norm(x-y);
		comp_N = @(x,y) -dot(x,y);
		% The image center position is the right choice to
		% 	discriminate as the wrong normal for two close isocontour
		% 	will also get close.
		[Xc_best,index_best] = get_best_ambiguous(Xc,nb_iso,2,comp_pos);
		N_best = get_best_ambiguous_index(N,nb_iso,3,index_best);
		N_best = N_best/norm(N_best);
	else
		Xc_best = Xc;
		N_best = N;
	end
	% Apply the convention for the normals
	N_best = convention_normals(N_best,Xc_best,K_im,pt_vis);
end
