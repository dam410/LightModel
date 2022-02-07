function [dssp1,dssp2,psp1,psp2,ps] = case_8_2p(data)

	% Known parameters

	% Estimation			
	[Xc,N] = plane_orientation_from_circular_contours(data.K,data.T_cam,...
		data.isocontour.CurveParameters{1});
	P_s = source_plane_from_circle_center_orientation(data.K,Xc,N);
	% TODO
	dssp = 0;
	psp = [N;0];
	ps = P_s;
end
