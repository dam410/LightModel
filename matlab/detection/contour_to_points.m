function [points,levels] = contour_to_points(C)
	points = {};
	levels = {};
	i_headreader = 1;
	i_current_contour = 0;
	n_data_contour = size(C,2);
	level_prev = 0;
	while i_headreader < n_data_contour
		nb_pts = C(2,i_headreader);
		level = C(1,i_headreader);
		if (level ~= level_prev)
			i_current_contour = i_current_contour+1;
			levels = [levels,level];
			points{i_current_contour} = C(:,(i_headreader+1):(i_headreader+nb_pts));
			level_prev = level;
		else
			points{i_current_contour} = [points{i_current_contour},...
				C(:,(i_headreader+1):(i_headreader+nb_pts))];
		end
		i_headreader = i_headreader+1+nb_pts;
	end
end
