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
		[Xc,N] = plane_orientation_from_circular_contours_co(data.K,data.T_cam,curves,pt_vis);
		%[Xc,N] = plane_orientation_from_circular_contours(data.K,data.T_cam,curves)
		% Distance to plane cannot be estimated, we set it to 1.0 and calculate h accordingly.
		psp{i_p} = {[N;0.0]};
		dssp{i_p} = {Inf};
	end
end
