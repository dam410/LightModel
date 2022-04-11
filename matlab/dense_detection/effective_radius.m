% Calculate the minimum and maximum radius associates with visible isocontours
%       points on the image
function [r_circle] = effective_radius(pt_in,H_opt,N_CIRC)
	% Test if the center is inside the polygon (then include it)
	P_center = inv(H_opt)*[0;0;1];
	p_center = P_center(1:2)./P_center(3);
	if inpolygon(p_center(1),p_center(2),pt_in(:,1),pt_in(:,2))
		pt_in = [pt_in;transpose(p_center)];
	end
        nb_pt = size(pt_in,1);
        proj_pt = [pt_in,ones(nb_pt,1)]*transpose(H_opt);
        proj_pt = proj_pt./proj_pt(:,3);
        r_pt = sqrt(proj_pt(:,1).^2+proj_pt(:,2).^2);
	% Prendre la médiane si on ne s'intéresse qu'à un seul isocontour
	if N_CIRC == 1
		r_circle = median(r_pt);
	else
		min_r = min(r_pt);
		max_r = max(r_pt);
		max_r = min_r+0.75*(max_r-min_r);
		r_circle_full = min_r:(max_r-min_r)/(N_CIRC+1):max_r;
		r_circle = r_circle_full(2:(end-1));
	end
end
