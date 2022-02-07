function [err] = eval_pose_1p(dssp,psp,ps,data)
	err_h = abs(dssp-data.groundtruth.ScenePlaneDistanceSource{1});
	N = psp(1:3);
	err_orient = 180*acos(dot(N,data.groundtruth.ScenePlaneOrientation{1}(:,3)))/pi;
	err_d = abs(psp(4)-data.groundtruth.ScenePlanePosition{1});
	% The evaluation of the source error depends on what we calculated
	if size(ps,1) == 3 && size(ps,2) == 1
		% We have a 3D point
		err_s = norm(ps-data.groundtruth.SourcePosition);
	elseif size(ps,1) == 4 && size(ps,2) == 1
		% We have a plane
		% Check normalization of the plane
		ps = ps/norm(ps(1:3));
		% Distance of the point to the plane
		err_s = dot(ps(1:3),data.groundtruth.SourcePosition)+ps(4);
	elseif size(ps,1)==4 && size(ps,2) == 4
		% We have a line in plucker matrix coordinates
		[d,m] = plucker_matrix_to_dm(ps);
		% Using the reciprocal product
		err_s = norm(cross_antisym(d)*data.groundtruth.SourcePosition+m);
	else
		disp('Not recognized source output');
		err_s = 0;
	end
	err = [err_h,err_orient,err_d,err_s];

	%TODO Deal with 2 planes

end
