function [err] = eval_pose_2(dssp1,dssp2,psp1,psp2,ps,data)
	err_h1 = abs(dssp1-data.groundtruth.ScenePlaneDistanceSource{1});
	err_h2 = abs(dssp2-data.groundtruth.ScenePlaneDistanceSource{2});
	N1 = psp1(1:3);
	N2 = psp2(1:3);
	err_orient1 = 180*acos(dot(N1,data.groundtruth.ScenePlaneOrientation{1}(:,3)))/pi;
	err_orient2 = 180*acos(dot(N2,data.groundtruth.ScenePlaneOrientation{2}(:,3)))/pi;
	err_d1 = abs(psp1(4)-data.groundtruth.ScenePlanePosition{1});
	err_d2 = abs(psp2(4)-data.groundtruth.ScenePlanePosition{2});
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
	err = [err_h1,err_h2,err_orient1,err_orient2,err_d1,err_d2,err_s];
end
