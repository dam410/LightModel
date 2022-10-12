% Special case 7 where the PLS is colocalised with camera
function [dssp,psp,ps] = case_7_co_ps(data)

	% Known parameters
	ps = {[0;0;0]};
	K = data.K;
	% Estimation
	if isfield(data,'isocontour')
		n_p = length(data.isocontour.CurveParameters);
	elseif isfield(data,'homography')
		n_p = length(data.isocontour.homography);
	else
		n_p = 1;
	end
	dssp = cell(1,n_p);
	psp = cell(1,n_p);
	for i_p = 1:n_p
		pt_vis = get_visible_point_from_data(data,i_p);
		if iscell(data.isocontour.CurveParameters{i_p}{1})
			curves = data.isocontour.CurveParameters{i_p}{1};
		else
			curves = data.isocontour.CurveParameters{i_p};
		end
		% Check if we actually only need to use the homography, because of top-down method
		skip = false;
		if isfield(data,'homography')
			if iscell(data.homography{i_p})
				curves = data.isocontour.CurveParameters{i_p}{1};
				H = inv(data.K)*inv(data.homography{i_p}{1});
			else
				curves = data.isocontour.CurveParameters{i_p};
				H = inv(data.K)*inv(data.homography{i_p});
			end
			% If H is only a rotation + canonical projection, we can use it directly to get N
			R = [H(:,1),H(:,2),H(:,3)/norm(H(:,3))];
			if norm(transpose(R)*R-eye(3))<1e-3
				N = R(:,3);
				Xc = N/N(3);
				skip = true;
			end
		if ~skip
			[Xc,N] = plane_orientation_from_circular_contours_co(data.K,data.T_cam,curves,pt_vis);
		end
		% Get two possible normal using a variation to plane_orientation_from_circular_contours_co
		%[Xc,N] = plane_orientation_from_circular_contours(data.K,data.T_cam,curves)
		% Distance to plane cannot be estimated, we set it to 1.0 and calculate h accordingly.
		psp{i_p} = {[N;0.0]};
		dssp{i_p} = {Inf};
	end
end
