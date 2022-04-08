% Update : Add an option for optimizing the same spline parameters for each plane
function [err_global,J_err_global] = evaluate_error_global_dense_derivate(...
		inv_K,...
		pts_cell,...
		I_pt_cell,...
		I_vect_cell,...
		param_vect,...
		display,...
		single_spline)
        S = param_vect(1:3);
        nb_plane = length(pts_cell);
        NB_R = length(I_vect_cell{1});
	nb_pts_total = sum(cellfun(@(x) length(x),pts_cell));
	if nargin < 7 || ~single_spline
		single_spline = false;
		i_param = 4;
	else
		i_param = 4+NB_R;
	end
	if nargout > 1
		if ~single_spline
			J_err_global = sparse(nb_pts_total,3+nb_plane*(3+NB_R));
		else
			J_err_global = sparse(nb_pts_total,3+nb_plane*(3));
		end
	end
	nb_pt_used = 1;
        err_global = zeros(nb_pts_total,1);
        for i_p = 1:nb_plane
		pts = pts_cell{i_p};
		n_pt_plane = length(pts);
		I_pt = I_pt_cell{i_p};
		% Depending if plane share the same spline, we extract parameters differently
		if ~single_spline
			I_vect = I_vect_cell{i_p};
			x_param = [S;param_vect(i_param:(i_param+2+NB_R))];

		else
			% If the single spline model is used, we only have one I_vect vector in I_vect_cell
			I_vect = I_vect_cell{1};
			x_param = [S;param_vect(i_param:(i_param+2));param_vect(4:3+NB_R)];
		end
			% NEW METHOD with Jacobian direct computation
		if nargout > 1 && single_spline
			[err,J_err] = evaluate_error_dense_plane_single_spline(inv_K,pts,...
				I_pt,I_vect,NB_R,x_param,display);
		elseif nargout > 1
			[err,J_err] = evaluate_error_dense_plane(inv_K,pts,I_pt,I_vect,NB_R,x_param,display);
		elseif single_spline
			[err] = evaluate_error_dense_plane_single_spline(inv_K,pts,...
				I_pt,I_vect,NB_R,x_param,display);
		else
			[err] = evaluate_error_dense_plane(inv_K,pts,I_pt,I_vect,NB_R,x_param,display);
		end
		% Then we rebuilt resulting Jacobian differently
		if ~single_spline && nargout > 1
			J_err_global(nb_pt_used:(nb_pt_used+n_pt_plane-1),[1:3,i_param:(i_param+2+NB_R)]) = J_err;
		elseif nargout>1
			J_err_global(nb_pt_used:(nb_pt_used+n_pt_plane-1),[1:3,i_param:(i_param+2)]) = J_err(:,1:6);
			J_err_global(nb_pt_used:(nb_pt_used+n_pt_plane-1),4:3+NB_R) = J_err(:,7:end);
		end
		if ~single_spline
			i_param = i_param+3+NB_R;
		else
			i_param = i_param+3;
		end
		err_global(nb_pt_used:(nb_pt_used+n_pt_plane-1)) = err;
		nb_pt_used = nb_pt_used + n_pt_plane;
        end
end

