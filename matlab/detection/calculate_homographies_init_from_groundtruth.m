function [homographies] = calculate_homographies_init_from_groundtruth(data)
	n_p = length(data.groundtruth.ScenePlaneOrientation);
	homographies = cell(1,n_p);
	for i_p = 1:n_p
		% Calculate the projection of the source on the plane
		X = data.groundtruth.SourcePosition+data.groundtruth.ScenePlaneDistanceSource{i_p}*data.groundtruth.ScenePlaneOrientation{i_p}(:,3);
		H = data.K*[data.groundtruth.ScenePlaneOrientation{i_p}(:,1:2),X];
		homographies{i_p} = H;
	end
end
