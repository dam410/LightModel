function [] = draw_line_image(p,q,size_I)
% Calculate the intersections points of the line with border
m = size_I(1);
n = size_I(2);
u = cross([p;1],[q;1]);
pt_corners = [1,1,n,n;1,m,m,1;1,1,1,1];
lines = cross(pt_corners(:,1:end),pt_corners(:,[2:end,1]));
pt_inters = cross([u,u,u,u],lines);
pt_inters = pt_inters(:,abs(pt_inters(3,:))>1e-3);
pt_inters = pt_inters./pt_inters(3,:);
dist_to_qp = dot(pt_inters(1:2,:)-p,repmat(q-p,1,size(pt_inters,2)));
ind_sup = find(dist_to_qp>0);
ind_inf = find(dist_to_qp<0);
[~,i_sup] = min(abs(dist_to_qp(ind_sup)));
[~,i_inf] = min(abs(dist_to_qp(ind_inf)));
pt_sup = pt_inters(1:2,ind_sup(i_sup));
pt_inf = pt_inters(1:2,ind_inf(i_inf));
line([pt_inf(1),pt_sup(1)],[pt_inf(2),pt_sup(2)],'Color','b');
end 

