function [table] = table_error_pose(err_h,err_orient,err_d,err_s)
	T.Error = {'Position (meter)';'Distance to source (meter)';'Orientation (Degree)'};
	T.Source = {num2str(err_s);'N/A';'N/A'};
	for i_p = 1:length(err_h)
		T.(genvarname(['Plane' num2str(i_p)])) = ...
			{num2str(err_d(i_p));num2str(err_h(i_p));num2str(err_orient(i_p))};
	end
	table = struct2table(T);
end
