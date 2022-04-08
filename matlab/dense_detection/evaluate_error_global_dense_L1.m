function [norm_L1,gradient_norm_L1,H_hessian_pattern] = evaluate_error_global_dense_L1(...
	inv_K,pts_cell,I_pt_cell,I_vect_cell,param_vect,display)

        S = param_vect(1:3);
        nb_plane = length(pts_cell);
        NB_R = length(I_vect_cell{1});
	nb_pts_total = sum(cellfun(@(x) length(x),pts_cell));
	% Initialize norm L1
	norm_L1 = 0;
	% Gradient vector
	gradient_norm_L1 = zeros(1,3+nb_plane*(3+NB_R));
	% Hessian pattern matrix
	H_hessian_pattern = sparse(3+nb_plane*(3+NB_R),3+nb_plane*(3+NB_R));
        i_param = 4;
	nb_pt_used = 1;
        for i_p = 1:nb_plane
		pts = pts_cell{i_p};
		n_pt_plane = length(pts);
		I_pt = I_pt_cell{i_p};
		I_vect = I_vect_cell{i_p};
		x_param = [S;param_vect(i_param:(i_param+2+NB_R))];
		% Calculate the local gradient for one plane
		[norm_L1_plane,gradient_norm_L1_plane] = evaluate_error_dense_plane_L1(...
			inv_K,pts,I_pt,I_vect,NB_R,x_param,display);
		norm_L1 = norm_L1 + norm_L1_plane;	
		gradient_norm_L1 = gradient_norm_L1 + gradient_norm_L1_plane;
		% Calculate the pattern of the Hessian
		H_hessian_pattern([1:3,i_param:(i_param+2+NB_R)],[1:3,i_param:(i_param+2+NB_R)]) = 1;
                if display==1
                        display_spline_function(r_vect,r_pt,I_vect_cell{i_p},I_pt_cell{i_p},['Curve radius/intensity i_p=',num2str(i_p)]);
                end
		i_param = i_param+3+NB_R;
		nb_pt_used = nb_pt_used + n_pt_plane;
        end
end


%%% OLD FUNCTIONS STUFF
		%save('data_test_evaluate_error_dense_plane.mat');
		%disp('data saved');
		%pause;
                %[N,J_N] = angle_to_normal(param_vect(i_param),param_vect(i_param+1));
		%% Calculate the Jacobian
		%J_N_full = [zeros(3),J_N,zeros(3,NB_R+1)];
		%J_S_full = [eye(3),zeros(3,NB_R+3)];
		%J_d = zeros(1,NB_R+6);
		%J_d(6) = 1;
                %d = param_vect(i_param+2);
                %r_vect = transpose(param_vect((i_param+3):(i_param+2+NB_R)));
                %i_param = i_param+3+NB_R;
                %nb_pt = length(pts_cell{i_p});
                %% Calculate the minimum and maximum radius using back projection on the plane
                %proj_pt = inv_K*transpose([pts_cell{i_p},ones(nb_pt,1)]);
                %% Calculate intersection with the plane
                %lambdas = -d./(transpose(N)*proj_pt);
                %proj_pt = lambdas.*proj_pt;
                %% Calculate projection of the source
                %h =  -d-dot(S,N);
		%J_h = -J_d-transpose(S)*J_N_full+transpose(N)*J_S_full;
                %Xs = S+h*N;
		%J_Xs = J_S_full + N.*J_h + h*J_N;
                %% Display the 3D points with the projection of the source
                %%figure('Name','3D points and source projection');
                %%plot3(Xs(1),Xs(2),Xs(3),'og');
                %%hold on;
                %%plot3(proj_pt(1,:),proj_pt(2,:),proj_pt(3,:),'+r');
                %%axis equal;
                %r_pt = transpose((proj_pt(1,:)-Xs(1)).^2+(proj_pt(2,:)-Xs(2)).^2+(proj_pt(3,:)-Xs(3)).^2);
                %I_pt_proj = spline(r_vect,I_vect_cell{i_p},r_pt);
