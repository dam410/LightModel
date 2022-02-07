% Initialize homography and radius control values for the spline
%	using points of maximal and minimal intensity
function [h_vect,r_vect] = initial_parameters_rough(I,pt_in_discrete,pt_in_poly,I_vect,I_pt)
	[n,m] = size(I);
	[I_max,ind_max] = max(I_pt);
        [I_min,ind_min] = min(I_pt);
        pt_max = [pt_in_poly(ind_max,1),pt_in_poly(ind_max,2)];
        pt_min = [pt_in_poly(ind_min,1),pt_in_poly(ind_min,2)];
        n_val = length(I_vect);
        r_max = norm(pt_min-pt_max)^2;
        r_min = 0;
	H = [1,0,-pt_max(1);0,1,-pt_max(2);0,0,1];
        h_vect = transpose(H(1:8));
	% Calculate using the projection of the points
	[pt_line_x,pt_line_y] = bresenham(pt_in_discrete(ind_max,2),pt_in_discrete(ind_max,1),...
		pt_in_discrete(ind_min,2),pt_in_discrete(ind_min,1));
	idx_pt_line = sub2ind([n,m],pt_line_x,pt_line_y);
	idx_pt_discrete = sub2ind([n,m],pt_in_discrete(:,2),pt_in_discrete(:,1));
	[~,ind_line] = ismember(idx_pt_line,idx_pt_discrete);
	ind_line = ind_line(ind_line>0);
	r_pt = radius_from_H(pt_in_poly,H);
	I_line = I_pt(ind_line);
	r_line = r_pt(ind_line);
	r_vect = [r_min];
	for i=2:(n_val-1)
		% Calculate the median but ensures the monotony
		ind_valid = find(I_line<=I_vect(i-1) & I_line>=I_vect(i+1) & r_line>r_vect(i-1));
		r_vect = [r_vect;mean(r_line(ind_valid))];
	end
	r_vect = [r_vect;r_max];
	% Problem avoid the non-monotonic values selection from r_pt
        %r_vect = transpose(r_min:(r_max-r_min)/(n_val-1):r_max);
end

