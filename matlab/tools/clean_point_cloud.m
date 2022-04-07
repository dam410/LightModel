% Remove the points to close to optical center manually
function [ptCloud_out] = clean_point_cloud(ptCloud,minDist)
	ind_select = find(sum(ptCloud.Location.^2,2) > minDist);
	Location = ptCloud.Location(ind_select,:);
	ptCloud_out = pointCloud(Location);
end
