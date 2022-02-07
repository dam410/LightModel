function [dssp,psp,ps] = visualize_results_surface_only(data,dssp,psp,ps,i_surface)

	if nargin < 5
		surf_index = 1:length(psp);
	else
		surf_index = [i_surface];
	end
	
	

	nb_scene_plane_estimated = length(psp);
	
	K_im = data.K;
	NB_PTS = 100;

	figure('Name','3D scene visualisation');
	hold on;

	% Show the camera
	R = [1 0 0; 0 1 0;0 0 1];
	O = zeros(3,1);
        t = transpose(O);
        pose = rigid3d(R,t);
        cam = plotCamera('AbsolutePose',pose,'Opacity',0.4,'Color','b','Size',0.05);


	% Show the surfaces and 100 random points of the isocontours for groundtruth and reconstruction
	for i_p=1:nb_scene_plane_estimated
		if isfield(data,'isocontour') && iscell(data.isocontour.Points{i_p}{1})
			n_ambig = length(data.isocontour.Points{i_p});
		else
			n_ambig = 1;
		end
		% Show the groundtruth surface planes but also 100 random points in the isocontours
		if isfield(data,'groundtruth') && isfield(data.groundtruth,'ScenePlanePolygon')
			points_3D = data.groundtruth.ScenePlanePolygon{surf_index(i_p)};
			points_3D = [points_3D,points_3D(:,1)];
			%patch(points_3D(1,:),points_3D(2,:),points_3D(3,:),[1,0,0],'FaceAlpha',.3,'FaceColor','m');
		else
			points_3D = inv(data.K)*[data.polys_2D{i_p};ones(1,size(data.polys_2D{i_p},2))];
		end
		% Display the reconstruction found
		for i_ambig = 1:length(psp{i_p})
			% Back project the points
			% Find the intersection of the points with the estimated plane
			lambdas = -psp{i_p}{i_ambig}(4)./(transpose(psp{i_p}{i_ambig}(1:3))*points_3D);
			points_3D_e = lambdas.*points_3D;
			points_3D_e = points_3D_e(:,[1:end,1]);
			%plot3(points_3D_e(1,:),points_3D_e(2,:),points_3D_e(3,:),'-b');
			patch(points_3D_e(1,:),points_3D_e(2,:),points_3D_e(3,:),[0,0,1],'FaceAlpha',.3);
			% If h=0 we set the observed point to a distance of 1 unit
			if norm(dssp{i_p}{i_ambig})<1e-6 && isfield(data,'isocontour') && isfield(data.isocontour,'Points')
				pt_avg = mean(data.isocontour.Points{i_p}{i_ambig}{1});
				ray = inv(data.K)*[pt_avg(2);pt_avg(1);1];
				ray = ray/norm(ray);
				psp{i_p}{i_ambig}(4) = -transpose(psp{i_p}{i_ambig}(1:3))*ray;
				dssp{i_p}{i_ambig} = distance_source_plane(psp{i_p}{i_ambig}(1:3),ps{1},psp{i_p}{i_ambig}(4));
			end
			% If homography is given for each plane in data, we use it to
			% calculate the projection of the source on the plane
			% and display the line coming from it to the source
			if isfield(data,'homography')
				H = data.homography{i_p}{i_ambig};
				invH = inv(H);
				Xc_3D = inv(data.K)*invH*[0;0;1];
				N = transpose(data.K)*transpose(H)*[0;0;1];
				N = N/norm(N);
				% Calculate its depth using the plane equation
				lambda = ...
					-psp{i_p}{i_ambig}(4)./(transpose(psp{i_p}{i_ambig}(1:3))*Xc_3D);
				Xc_3D = lambda*Xc_3D;
			else
				Xc_3D = ps{1}+dssp{i_p}{i_ambig}*psp{i_p}{i_ambig}(1:3);
				%Xc_end = Xc_3D-0.9*dssp{i_p}{i_ambig}*psp{i_p}{i_ambig}(1:3);
			end
			Xc_end = -0.92*dssp{i_p}{i_ambig}*psp{i_p}{i_ambig}(1:3);
			%quiver3(Xc_3D(1),Xc_3D(2),Xc_3D(3),Xc_end(1),Xc_end(2),Xc_end(3),'MarkerMode','manual','Color','r',...
			%	'MaxHeadSize',0.2,'Linewidth',2,'MarkerSize',5,'Autoscale','off','AutoScaleFactor',2);
		end
	end

	%% Show the source
	%if isfield(data,'groundtruth')
	%	S = data.groundtruth.SourcePosition;
	%	plot3(S(1),S(2),S(3),'pk','MarkerFaceColor','m','MarkerSize',15);
	%end
	%for i_ambig = 1:length(ps)
	%	S_est = ps{i_ambig};
	%	if size(S_est,1) == 3 && size(S_est,2) == 1
	%		% We have a 3D point
	%		plot3(S_est(1),S_est(2),S_est(3),'pk','MarkerFaceColor','y','MarkerSize',15);
	%	elseif size(S_est,1) == 4 && size(S_est,2) == 1
	%		% We have a plane
	%		% We can display the two orthogonal axis passing through the camera center
	%		UV = null(S_est(1:3)*transpose(S_est))
	%		plot3([0,UV(1,1)],[0,UV(2,1)],[0,UV(3,1)],'-r');
	%		plot3([0,UV(1,2)],[0,UV(2,2)],[0,UV(3,2)],'-r');
	%	elseif size(S_est,1)==4 && size(S_est,2) == 4
	%		% We have a line in plucker matrix coordinates
	%		[d,m] = plucker_matrix_to_dm(S_est);
	%		m = m/norm(d);
	%		d = d/norm(d);
	%		% Find the closest point on the line to the source
	%		m_S = -m - cross(S,d);
	%		X_S = S + cross(d,m_S);

	%		plot3([X_S(1)+3*d(1),X_S(1)-3*d(1)],...
	%			[X_S(2)+3*d(2),X_S(2)-3*d(2)],...
	%			[X_S(3)+3*d(3),X_S(3)-3*d(3)],'-r');
	%	else
	%		disp('Not recognized source output');
	%	end
	%end
	axis equal;
	set(gca,'visible','off')
end

