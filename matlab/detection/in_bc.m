% This function take two points and the frame of an image that is a boundary
%	* It returns true if the points is inside the boundary
%	* Otherwise it returns the intersection with the boundary
function [inboundcheck,pt_inter] = in_bc(pt,size_I,pt_from)
        n = size_I(1);
        m = size_I(2);
        inboundcheck = (pt(1) > 0) && (pt(1)<(n+1)) && (pt(2) > 0) && (pt(2) < (m+1));
        if nargin==3
                if ~inboundcheck
                        seg_up = [1,1;1,m];
                        seg_down = [n,n;1,m];
                        seg_right = [1,n;m,m];
                        seg_left = [1,n;1,1];
			pt_inter_up = seg_inter([pt_from,pt],seg_up);
			pt_inter_down = seg_inter([pt_from,pt],seg_down);
			pt_inter_left = seg_inter([pt_from,pt],seg_left);
			pt_inter_right = seg_inter([pt_from,pt],seg_right);
                        pts_inter = [pt_inter_up,pt_inter_down,pt_inter_left,pt_inter_right];
                        % Take the closest point to pt
                        diff = pt_from-pts_inter;
                        dist_2 = diff(1,:).^2+diff(2,:).^2;
                        [~,ind_closest] = min(dist_2);
                        pt_inter = pts_inter(:,ind_closest);
                else
                        pt_inter = [];
                end
end

