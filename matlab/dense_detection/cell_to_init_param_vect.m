% This function is the preliminary step before global optimization
%	* Calculate the initial parameters vector with our data
%	* Calculate the Jacobian pattern matrix for least square problem
%	* Calculate the Hessian pattern matrix for the L1 or L2 problem
%	* Calculate the constraint vector to avoid non monotonic spline behavior
% Update : Add an option for optimizing the same spline parameters for each plane
function [param_vect,J_sparse,H_sparse,constraints_vector] = cell_to_init_param_vect(inv_K,ps,psp,pts_cell,NB_R,single_spline)
	if nargin < 6
		single_spline = false;
	end
	S = ps{1};
	EPS = 1;
        nb_plane = length(psp);
        nbs_pt = ones(1,nb_plane);
	R_vect = tril(ones(NB_R));
        inv_R_vect = inv(R_vect);
	constraints_vector = -Inf*ones(3,1);
	if ~single_spline
		param_vect = zeros(3+(3+NB_R)*nb_plane,1);
	else
		param_vect = zeros(3+NB_R+3*nb_plane,1);
		MIN_rpt = Inf;
		MAX_rpt = -Inf;
		constraints_vector = [constraints_vector(:);-Inf;zeros(NB_R-1,1)];
	end
	param_vect(1:3) = S;
        for i_p = 1:nb_plane
                nbs_pt(i_p) = length(pts_cell{i_p});
		% If no initial psp is found, we use classical plane pose d = -1 and N = [0;0;1]
		if isempty(psp{i_p})
			N = [0;0;1];
			d = -1;
		else
			N = psp{i_p}{1}(1:3);
			d = psp{i_p}{1}(4);
		end
		[alpha,beta] = normal_to_angle(N);
                % Calculate the minimum and maximum radius using back projection on the plane
                proj_pt = inv_K*transpose([pts_cell{i_p},ones(nbs_pt(i_p),1)]);
                %proj_pt = proj_pt./proj_pt(3,:);
                % Calculate intersection with the plane
                lambdas = -d./(transpose(N)*proj_pt);
                proj_pt = lambdas.*proj_pt;
		% If some assumption can be made on the attenuation function, we interpolate a new function
		if single_spline
			% Calculate their distance to S
			r_S = sqrt((proj_pt(1,:)-S(1)).^2 + (proj_pt(2,:)-S(2)).^2 + (proj_pt(3,:)-S(3)).^2);
			% Calculate the angle between direction to S and N for each point
			cos_thetas = dot((proj_pt-S)./r_S,N*ones(1,nbs_pt(i_p)));
			% Calculate the new r_pt using our prior on attenuation function
			r_pt = cos_thetas./(r_S.^2+ EPS);
		else
			% Calculate projection of the source
			h =  -d-dot(S,N);
			Xs = S+h*N;
			r_pt = transpose((proj_pt(1,:)-Xs(1)).^2+(proj_pt(2,:)-Xs(2)).^2+(proj_pt(3,:)-Xs(3)).^2);
		end
                max_rpt = max(r_pt);
                min_rpt = min(r_pt);
		if ~single_spline
			r_vect = transpose(min_rpt:(max_rpt-min_rpt)/(NB_R-1):max_rpt);
			a_vect = inv_R_vect*r_vect;
			param_vect(4+(i_p-1)*(3+NB_R):3+(i_p)*(3+NB_R)) = [alpha;beta;d;a_vect];
			% Update the constraints vector
			constraints_vector = [constraints_vector(:);-Inf*ones(4,1);zeros(NB_R-1,1)];
		else
			param_vect(4+NB_R+(i_p-1)*3:3+NB_R+(i_p)*3) = [alpha;beta;d];
			% Update the constraints vector
			constraints_vector = [constraints_vector(:);-Inf*ones(3,1)];
			MIN_rpt = min(min_rpt,MIN_rpt);
			MAX_rpt = max(max_rpt,MAX_rpt);
		end
        end
	if ~single_spline
		J_sparse = sparse(sum(nbs_pt),3+(3+NB_R)*nb_plane);
		% Jacobian columns for the source
		J_sparse(:,1:3) = 1;
		nb_pt = 1;
		for i_p = 1:nb_plane
			J_sparse(nb_pt:(nb_pt-1+nbs_pt(i_p)),(4+(i_p-1)*(3+NB_R)):(3+(i_p)*(3+NB_R))) = 1;
			%H_sparse([1:3,(4+(i_p-1)*(3+NB_R)):(3+(i_p)*(3+NB_R))],[1:3,(4+(i_p-1)*(3+NB_R)):(3+(i_p)*(3+NB_R))]) = 1;
			nb_pt = nb_pt+nbs_pt(i_p);
		end
		H_sparse = double((transpose(J_sparse)*J_sparse>0));
	else
		param_vect(4:(3+NB_R)) = inv_R_vect*transpose(MIN_rpt:(MAX_rpt-MIN_rpt)/(NB_R-1):MAX_rpt);
		J_sparse = sparse(sum(nbs_pt),3+NB_R+3*nb_plane);
		% Jacobian columns for the source	
		J_sparse(:,1:3+NB_R) = 1;
		nb_pt = 1;
		for i_p = 1:nb_plane
			J_sparse(nb_pt:(nb_pt-1+nbs_pt(i_p)),4+NB_R+(i_p-1)*3:(3+NB_R+i_p*3)) = 1;
			nb_pt = nb_pt+nbs_pt(i_p);
		end
		H_sparse = double((transpose(J_sparse)*J_sparse>0));
	end
end

