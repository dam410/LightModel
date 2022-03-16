function [homographies,isocontours] = get_groundtruth_isocontours(data,polys_2D)

	n_p = length(data.groundtruth.ScenePlaneOrientation);
	homographies = cell(1,n_p);
	isocontours = struct();
	isocontours.Points = cell(1,n_p);
	isocontours.CurveParameters = cell(1,n_p);
	for i_p=1:n_p
		X = data.groundtruth.SourcePosition - ...
			data.groundtruth.ScenePlaneDistanceSource{i_p}*...
			data.groundtruth.ScenePlaneOrientation{i_p}(:,3);
		H = data.K*[data.groundtruth.ScenePlaneOrientation{i_p}(:,1:2),X];
		invH = inv(H);
		homographies{i_p} = invH;
		if isfield(data,'isocontour')
			n_iso = length(data.isocontour.Points{i_p}{1});
			[points_iso,curves_iso] = isocontours_from_homography(polys_2D{i_p},invH,n_iso);
			points_iso
			isocontours.Points{i_p} = points_iso;
			isocontours.CurveParameters{i_p} = curves_iso;
		end
	end
	
end
