% This function find linear separator in angle space to keep two solution for the normal appart
function [a_1,a_2,lb_1,ub_1,lb_2,ub_2] = normal_angle_separation(a_1,b_1,a_2,b_2)
	lb_1 = [-Inf;0];
	ub_1 = [Inf;pi/2];
	lb_2 = [-Inf;0];
	ub_2 = [Inf;pi/2];
	% Only put a separator on b angle if the differenc
	a_21 = atan2(sin(a_2-a_1),cos(a_2-a_1));
	b_21 = atan2(sin(b_2-b_1),cos(b_2-b_1));
	perc_sep = 0.4;
	if abs(a_21)>pi/8
		a_2 = a_1+a_21;
		lb_1(1) = a_1-perc_sep*abs(a_21);
		ub_1(1) = a_1+perc_sep*abs(a_21);
		lb_2(1) = a_2-perc_sep*abs(a_21);
		ub_2(1) = a_2+perc_sep*abs(a_21);
	elseif abs(b_21)>pi/8
		lb_1(2) = b_1-perc_sep*abs(b_21);
		ub_1(2) = b_1+perc_sep*abs(b_21);
		lb_2(2) = b_2-perc_sep*abs(b_21);
		ub_2(2) = b_2+perc_sep*abs(b_21);
	else
		disp('Les normales sont trop proches pour être séparés');
	end

end
