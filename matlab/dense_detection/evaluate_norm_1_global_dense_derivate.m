function [norm_L1,gradient_norm_L1] = evaluate_norm_1_global_dense_derivate(inv_K,pts_cell,I_pt_cell,I_vect_cell,x_param,single_spline)
        if nargin < 6
                single_spline = false;
        end
        if nargout > 1
                [err,J_err] = evaluate_error_global_dense_derivate(inv_K,pts_cell,I_pt_cell,I_vect_cell,x_param,0,single_spline);
                [norm_L1,gradient_norm_L1] = error_norm_L1(err,J_err);
        else
                err = evaluate_error_global_dense_derivate(inv_K,pts_cell,I_pt_cell,I_vect_cell,x_param,0,single_spline);
                [norm_L1] = error_norm_L1(err);
        end
end
