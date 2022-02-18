% Extract a visible point for each scene plane
function [pts_vis] = get_visible_point_from_data(data,i_p)
	if isfield(data,'isocontour') && isfield(data.isocontour,'Points')
		if iscell(data.isocontour.Points{i_p}{1})
			pts_vis = data.isocontour.Points{i_p}{1}{1}(1,:);
		else
			pts_vis = data.isocontour.Points{i_p}{1}(1,:);
		end
	else
		pts_vis = [];
	end
end
