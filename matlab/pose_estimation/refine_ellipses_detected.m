function [curve_params_out] = refine_ellipses_detected_affine(isocontour_pts,curve_params)
	% Adding constraints to the optimization: 
	%	* their centers are closed to each others
	%	* No intersection
	MAX_PTS = 100;
	n_p = length(isocontour_pts);
	n_iso = length(isocontour_pts{1});
	curve_params_out = cell(1,n_p);
	for i_p = 1:2
		% Create an optimization vectors with ellipse's parameters
		params_ell_vect = cell2mat(curve_params{i_p});
		angles = [];
		pts_kept = cell(1,n_iso);
		for i_iso = 1:n_iso
			nb_pts = length(isocontour_pts{i_p}{i_iso});	
			if nb_pts>MAX_PTS
				samples = datasample(1:nb_pts,MAX_PTS,'Replace',false);
			else
				samples = datasample(1:nb_pts,MAX_PTS,'Replace',true);
			end
			pts_kept{i_iso} = isocontour_pts{i_p}{i_iso}(samples,:);
			pts_centered = pts_kept{i_iso} - curve_params{i_p}{i_iso}(1:2);
			cos_t = cos(curve_params{i_p}{i_iso}(5));
			sin_t = sin(curve_params{i_p}{i_iso}(5));
			a = real(curve_params{i_p}{i_iso}(3));
			b = real(curve_params{i_p}{i_iso}(4));
			x_0 = curve_params{i_p}{i_iso}(1);
			y_0 = curve_params{i_p}{i_iso}(2);
			xpts = (cos_t.*pts_centered(:,1)+...
				sin_t.*pts_centered(:,2))./a;
			ypts = (-sin_t.*pts_centered(:,1)+...
				cos_t.*pts_centered(:,2))./b;
			angle_loc = atan2(ypts,xpts);

			%t = 0:0.05:(2*pi);
			%xpts_ell = cos_t.*(a*cos(t))-sin_t.*(b*sin(t))+curve_params{i_p}{i_iso}(1);
			%ypts_ell = sin_t.*(a*cos(t))+cos_t.*(b*sin(t))+curve_params{i_p}{i_iso}(2);

			%xpts_ell = cos_t.*(a*cos(angle_loc))-sin_t.*(b*sin(angle_loc))+x_0;
			%ypts_ell = sin_t.*(a*cos(angle_loc))+cos_t.*(b*sin(angle_loc))+y_0;
			%figure('Name','Plot the ellipse with many points');
			%plot(xpts_ell,ypts_ell,'+b');
			%hold on;
			%plot(pts_kept{i_iso}(:,1),pts_kept{i_iso}(:,2),'+r');
			%figure('Name','Rectify on circles');
			%plot(xpts,ypts,'+r');
			%hold on;
			%plot(cos(angle_loc),sin(angle_loc),'+b');
			%axis equal;
			%figure('Name','Look points on the circle');
			%plot(xpts,ypts,'+k');
			%axis equal;
			%curve_params{i_p}{i_iso}(3)
			%curve_params{i_p}{i_iso}(4)
			angles = [angles,transpose(angle_loc)];
		end
		% Optimizing
		x_init = [params_ell_vect,angles];
		err = evaluate_ellipses(x_init,n_iso,pts_kept);
		fun_opt = @(param) evaluate_ellipses(param,n_iso,pts_kept);
		options = optimoptions('lsqnonlin');
		options.Algorithm = 'levenberg-marquardt';
                options.Display = 'iter-detailed';
                options.UseParallel = true;
                options.MaxIterations = 30;
		tic;
		x_opt = lsqnonlin(fun_opt,x_init,[],[],options);
		toc
		param_ell_vect_out = reshape(x_opt(1:(5*n_iso)),5,n_iso);
		curve_params_out{i_p} = cell(1,n_iso);
		for i_iso = 1:n_iso
			curve_params_out{i_p}{i_iso} = transpose(param_ell_vect_out(:,i_iso));
		end
	end
end

% params_vect contains the ellipses parameters and the
%	angles of the corresponding closest points on the ellipse
function [err] = evaluate_ellipses(params_vect,n_iso,points)
	params_ell_vect = reshape(params_vect(1:(5*n_iso)),5,n_iso);
	i_current = 5*n_iso;
	err = [];
	for i_iso=1:n_iso
		nb_pt = length(points{i_iso});
		angles = transpose(params_vect((1+i_current):(i_current+nb_pt)));
		cos_t = cos(params_ell_vect(5,i_iso));
		sin_t = sin(params_ell_vect(5,i_iso));
		a = params_ell_vect(3,i_iso);
		b = params_ell_vect(4,i_iso);
		x_0 = params_ell_vect(1,i_iso);
		y_0 = params_ell_vect(2,i_iso);
		xpts_ell = cos_t.*(a*cos(angles))-sin_t.*(b*sin(angles))+x_0;
		ypts_ell = sin_t.*(a*cos(angles))+cos_t.*(b*sin(angles))+y_0;

		%xpts = params_ell_vect(3,i_iso)*cos(angles);
		%ypts = params_ell_vect(4,i_iso)*sin(angles);
		%Xpts = transpose(params_ell_vect(1:2,i_iso))+...
		%	[cos(params_ell_vect(5,i_iso))*xpts-sin(params_ell_vect(5,i_iso))*ypts,...
		%	sin(params_ell_vect(5,i_iso))*xpts+cos(params_ell_vect(5,i_iso))*ypts];
		Xpts = [xpts_ell,ypts_ell];
		diff_pts = points{i_iso} - Xpts;
		err = [err;diff_pts(:)];
		%figure('Name','Show the points and difference');
		%plot(transpose([params_ell_vect(1,i_iso)*ones(nb_pt,1),points{i_iso}(:,1)]),...
		%transpose([params_ell_vect(2,i_iso)*ones(nb_pt,1),points{i_iso}(:,2)]),'-k');
		%hold on;
		%plot(transpose([params_ell_vect(1,i_iso)*ones(nb_pt,1),Xpts(:,1)]),...
		%transpose([params_ell_vect(2,i_iso)*ones(nb_pt,1),Xpts(:,2)]),'-k');
		%plot(Xpts(:,1),Xpts(:,2),'+r');
		%plot(points{i_iso}(:,1),points{i_iso}(:,2),'+b');
		i_current = i_current+length(points{i_iso});
		% Calculate the distance between the current center and all other centers
	end
end
