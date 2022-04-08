function [J_err_global] = jacobian_pattern_global_dense(pts_cell,I_vect_cell)
        nb_plane = length(pts_cell);
        NB_R = length(I_vect_cell{1});
	nb_pts_total = sum(cellfun(@(x) length(x),pts_cell));
	J_err_global = sparse(nb_pts_total,3+nb_plane*(3+NB_R));
        i_param = 4;
	nb_pt_used = 1;
        for i_p = 1:nb_plane
		pts = pts_cell{i_p};
		n_pt_plane = length(pts);
		J_err_global(nb_pt_used:(nb_pt_used+n_pt_plane-1),[1:3,i_param:(i_param+2+NB_R)]) = 1;
		i_param = i_param+3+NB_R;
		nb_pt_used = nb_pt_used + n_pt_plane;
        end
end
