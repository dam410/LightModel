function [S,psp,r_vects] = param_vect_to_cell(param_vect,nb_plane,NB_R,single_spline)
        S = param_vect(1:3);
        psp = cell(1,nb_plane);
        r_vects = cell(1,nb_plane);
        J_r_vect = tril(ones(NB_R));
        if nargin<4 || ~single_spline
                single_spline = false;
                i_param = 4;
        else
                i_param = 4+NB_R;
                a_vect = transpose(J_r_vect)*param_vect(4:(3+NB_R));
        end
        for i_p = 1:nb_plane
                if ~single_spline
                        N = angle_to_normal(param_vect(i_param),param_vect(i_param+1));
                        d = param_vect(i_param+2);
                        psp{i_p} = {[N;d]};
                        r_vects{i_p} = transpose(J_r_vect*param_vect((i_param+3):(i_param+2+NB_R)));
                        i_param = i_param+3+NB_R;
                else
                        N = angle_to_normal(param_vect(i_param),param_vect(i_param+1));
                        d = param_vect(i_param+2);
                        psp{i_p} = {[N;d]};
                        r_vects{i_p} = a_vect;
                        i_param = i_param+3;
                end
        end
end
