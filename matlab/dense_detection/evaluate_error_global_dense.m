function [err] = evaluate_error_global_dense(inv_K,pts_cell,I_pt_cell,I_vect_cell,param_vect,display)
	S = param_vect(1:3);
	nb_plane = length(pts_cell);
	NB_R = length(I_vect_cell{1});
	err = [];
	i_param = 4;
	for i_p = 1:nb_plane
		N = angle_to_normal(param_vect(i_param),param_vect(i_param+1));
		d = param_vect(i_param+2);
		r_vect = transpose(param_vect((i_param+3):(i_param+2+NB_R)));
		i_param = i_param+3+NB_R;
		nb_pt = length(pts_cell{i_p});
		% Calculate the minimum and maximum radius using back projection on the plane
		proj_pt = inv_K*transpose([pts_cell{i_p},ones(nb_pt,1)]);
		% Calculate intersection with the plane
		lambdas = -d./(transpose(N)*proj_pt);
		proj_pt = lambdas.*proj_pt;
		% Calculate projection of the source
		h =  -d-dot(S,N);
		Xs = S+h*N;
		% Display the 3D points with the projection of the source
		%figure('Name','3D points and source projection');
		%plot3(Xs(1),Xs(2),Xs(3),'og');
		%hold on;
		%plot3(proj_pt(1,:),proj_pt(2,:),proj_pt(3,:),'+r');
		%axis equal;
		r_pt = transpose((proj_pt(1,:)-Xs(1)).^2+(proj_pt(2,:)-Xs(2)).^2+(proj_pt(3,:)-Xs(3)).^2);
		I_pt_proj = spline(r_vect,I_vect_cell{i_p},r_pt);
		err = [err;(I_pt_proj-I_pt_cell{i_p})];
		if display==1
			display_spline_function(r_vect,r_pt,I_vect_cell{i_p},...
				I_pt_cell{i_p},['Curve radius/intensity i_p=',num2str(i_p)]);
		end
	end
end

