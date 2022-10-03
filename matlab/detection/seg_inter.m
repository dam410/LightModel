% Calculate the intersection between two segments return an empty vector if there is none
%	* Segments are given with their two extremity points
%	* s = [x0,x1;y0,y1]
function [pt] = seg_inter(seg_1,seg_2)
	EPS = 1e-10;
	p = seg_1(:,1);
	q = seg_2(:,1);
	r = seg_1(:,2)-p;
	s = seg_2(:,2)-q;
	% Calculate using  information given in the following link
	% https://stackoverflow.com/questions/563198/how-do-you-detect-where-two-line-segments-intersect
	cross_2D_rs = cross_2D(r,s);
	cross_2D_p_qr = cross_2D(p-q,r);
	if abs(cross_2D_rs)<EPS && abs(cross_2D_p_qr)<EPS
		t_0 = dot(q-p,r)/dot(r,r);
		t_1 = dot(q+s-p,r)/dot(r,r);
		t = sort([t_0,t_1]);
		[sort_t,check_t] = sort([0,1,t(1),t(2)]);
		check_inter = (check_t(1)==1 & check_t(2)==2) | (check_t(1)==3 & check_t(2)==4);
		if check_t
			pt = [];
		else
			pt = p+sort_t(2)*r;
		end
	elseif abs(cross_2D_rs)<EPS
		pt = [];
	else 
		t = cross_2D(q-p,s)/cross_2D(r,s);
		u = cross_2D(q-p,r)/cross_2D(r,s);
		if in_01([t,u])
			pt = [p+t*r];
		else
			pt = [];
		end
	end
end

function [result] = in_01(T)
	check = T>=0 & T<=1;
	result = prod(check);
end

function [res] = cross_2D(p,q)
	res = p(1)*q(2)-p(2)*q(1);
end
