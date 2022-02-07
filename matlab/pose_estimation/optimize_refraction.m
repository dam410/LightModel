% This function compute the intersection of the point with the miror using the refraction model
function [X,n_1,h,err] = optimize_refraction(inters,normals)
	norm_inters = sqrt(sum(inters.^2));
	inters_dir = inters./norm_inters;
	dir_ortho = (inters-dot(inters,normals));
	norm_dir_ortho = sqrt(sum(dir_ortho.^2));
	dir_ortho_normalized = dir_ortho./norm_dir_ortho;
	cos_i0 = dot(inters_dir,normals);
	sin_i0 = sqrt(1-cos_i0.^2);
	tan_i0 = sin_i0./cos_i0;
	sin_i1 = 1/r*sin_i0;
	tan_i1 = 1./(sqrt(1-sin_i1.^2));
	a_i = h*tan_i1;
	b_i = a_i./tan_i0;
	t_i = a_i.*dir_ortho_normalized+b_i.*normals;
	dirs = 2*inters-2*dot(inters,normals(1:3,:)).*normals(1:3,:);
	x_i1 = inters+t_i;
	x_i2 = dirs+t_i;
	figure('Name','Source line direction and intersection point');
	plot3([inters(1,:);dirs(1,:)],...
		[inters(2,:);dirs(2,:)],...
		[inters(3,:);dirs(3,:)],...
		'-r');
	hold on;
	plot3([x_i1(1,:);x_i2(1,:)],...
		[x_i1(2,:);x_i2(2,:)],...
		[x_i1(3,:);x_i2(3,:)],...
		'-b');
	L_cell = lines_2p_to_plucker(inters,dirs);
	[X_1,err] = closest_point_lines(L_cell)

	L_cell = lines_2p_to_plucker(x_i1,x_i2);
	[X_2,err] = closest_point_lines(L_cell)
	plot3(X_1(1),X_1(2),X_1(3),'+r','MarkerSize',10);
	plot3(X_2(1),X_2(2),X_2(3),'+b','MarkerSize',10);
end

function [err] = calculate_error(inters,normals,h,r)
        norm_inters = sqrt(sum(inters.^2));
        inters_dir = inters./norm_inters;
        dir_ortho = (inters-dot(inters,normals));
        norm_dir_ortho = sqrt(sum(dir_ortho.^2));
        dir_ortho_normalized = dir_ortho./norm_dir_ortho;
        cos_i0 = dot(inters_dir,normals);
        sin_i0 = sqrt(1-cos_i0.^2);
        tan_i0 = sin_i0./cos_i0;
        sin_i1 = 1/r*sin_i0;
        tan_i1 = 1./(sqrt(1-sin_i1.^2));
        a_i = h*tan_i1;
        b_i = a_i./tan_i0;
        t_i = a_i.*dir_ortho_normalized+b_i.*normals;
        dirs = 2*inters-2*dot(inters,normals(1:3,:)).*normals(1:3,:);
        x_i1 = inters+t_i;
        x_i2 = dirs+t_i;
        L_cell = lines_2p_to_plucker(x_i1,x_i2);
        [X,err] = closest_point_lines(L_cell)
end
